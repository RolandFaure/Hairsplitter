import pandas as pd 
import numpy as np
import sys

import gurobipy as grb
import pysam as ps

from sklearn.cluster import FeatureAgglomeration
from sklearn.cluster import AgglomerativeClustering
from sklearn.impute import KNNImputer
from sklearn.metrics import pairwise_distances

import time
from argparse import ArgumentParser

def get_data(file, contig_name,start_pos,stop_pos):
    #INPUT: a SORTED and INDEXED Bam file
    # Go through a window w on a contig and select suspicious positions
    itercol = file.pileup(contig = contig_name,start = start_pos,stop = stop_pos,truncate=True,min_base_quality=10)
    list_of_sus_pos = {}

    for pileupcolumn in itercol:
        if pileupcolumn.nsegments >=5:
            
            tmp_dict = {}
            sequence = np.char.upper(np.array(pileupcolumn.get_query_sequences()))
            bases, freq = np.unique(np.array(sequence),return_counts=True)
    
            if bases[0] =='':    
                ratio = freq[1:]/sum(freq[1:])
                bases = bases[1:]
            else:
                ratio = freq/sum(freq)
            
            if len(ratio) >0:
                idx_sort = np.argsort(ratio)
                
                if ratio[idx_sort[-1]]<0.95:
                    for pileupread in pileupcolumn.pileups:        
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.
                            if pileupread.alignment.query_sequence[pileupread.query_position] == bases[idx_sort[-1]]:
                                tmp_dict[pileupread.alignment.query_name] =  1
                            elif pileupread.alignment.query_sequence[pileupread.query_position] == bases[idx_sort[-2]]:
                                tmp_dict[pileupread.alignment.query_name] =  0
                            else :
                                tmp_dict[pileupread.alignment.query_name] = np.nan
                
                    list_of_sus_pos[pileupcolumn.reference_pos] = tmp_dict        
    return list_of_sus_pos

def pre_processing(X_matrix, min_col_quality = 3):
    ###Filling the missing values using KNN
    m,n = X_matrix.shape
    if m>1 and n>1:
        imputer = KNNImputer(n_neighbors= 10)
        matrix = imputer.fit_transform(X_matrix)
        upper,lower = 0.7,0.3
        matrix[(matrix>=upper)] = 1
        matrix[(matrix<=lower)] = 0
        matrix[(matrix<upper) * (matrix>lower)] = -1
    else:
        matrix = X_matrix  
    
    ###Split in to distinct regions(columns)
    regions = [list(range(matrix.shape[1]))]
    inhomogenious_regions = [list(range(matrix.shape[1]))]
    steps = []
    
    if m>5 and n > 15:
        agglo=FeatureAgglomeration(n_clusters=None,metric = 'hamming', linkage = 'complete',distance_threshold=0.35)
        agglo.fit(matrix)
        labels = (agglo.labels_)
        splitted_cols = {}
    
        for idx,label in enumerate(labels):
            if label in splitted_cols:
                splitted_cols[label].append(idx)
            else:
                splitted_cols[label] = [idx]
        
        regions = []        
        for cols in splitted_cols.values():
            
            ##if the region is too large, remove the noisiest reads by splitting them with a strict distance
            ##threshold, then take only the significant clusters
            if len(cols) > 15:
                matrix_reg = matrix[:,cols].copy()
                agglo=FeatureAgglomeration(n_clusters=None,metric = 'hamming', linkage = 'complete',distance_threshold=0.025)
                agglo.fit(matrix_reg)
                labels = (agglo.labels_)
                groups = {}

                for idx,label in enumerate(labels):
                    if label in groups:
                        groups[label].append(cols[idx])
                    else:
                        groups[label] = [cols[idx]]

                for c in groups.values():
                    if len(c) >= min_col_quality:
                        regions.append(c)
            elif len(cols)>=min_col_quality:
                regions.append(cols)
                
        inhomogenious_regions = []        
        for region in regions:
            thres = (min_col_quality/len(region))
            matrix_reg = matrix[:,region].copy()
            matrix_reg[matrix_reg==-1] = 0
            x_matrix = (matrix_reg.sum(axis = 1))/len(region)
            if len(x_matrix[(x_matrix>=thres)*(x_matrix<=1-thres)]) > 10:
                inhomogenious_regions.append(region)
            else:
                ###cut into 2 regions of 1 and 0
                agglo = AgglomerativeClustering(n_clusters = 2, metric = 'hamming', linkage = 'complete')
                agglo.fit(matrix_reg)  
                labels = agglo.labels_  
                cluster1,cluster0 = [],[]
                for idx,label in enumerate(labels):
                    if label == 0:
                        cluster0.append(idx)
                    else:
                        cluster1.append(idx)
                mat_cl1 = matrix_reg[cluster1,:].sum()/(len(cluster1)*len(region))
                mat_cl0 = matrix_reg[cluster0,:].sum()/(len(cluster0)*len(region))

                if (mat_cl0 > 0.9 or mat_cl0<0.1) and (mat_cl1 > 0.9 or mat_cl1<0.1):
                    steps.append((cluster1,cluster0,region))
            
    return matrix, inhomogenious_regions, steps

def quasibiclique(X_matrix, error_rate = 0.025):
    #Finding quasibiclique of a binary matrix
    X_problem = X_matrix.copy()
    
    cols_sorted = np.argsort(X_problem.sum(axis = 0))[::-1]
    rows_sorted = np.argsort(X_problem.sum(axis = 1))[::-1]
    
    m = len(rows_sorted)
    n = len(cols_sorted)
    if m==0 or n==0:
        return ([],[],False)
    #SELECTING MOSTLY 1 REGION FOR SEEDING
    seed_rows = m//3
    seed_cols = n//3
    if n>50:
        step_n = 10
    else:
        step_n = 2
    for x in range(m//3,m,10):
        for y in range(n//3,n,step_n):
            nb_of_1_in_seed = 0
            for i in rows_sorted[:x]:
                for j in cols_sorted[:y]:
                    nb_of_1_in_seed = nb_of_1_in_seed + X_problem[i][j]
            ratio_of_1 = nb_of_1_in_seed/(x*y)
            
            if ratio_of_1 > 0.99 and x*y>seed_rows*seed_cols:
                seed_rows = x
                seed_cols = y
                
    model = grb.Model('max_model')
    model.Params.OutputFlag = 0
    model.Params.MIPGAP = 0.05
    model.Params.TimeLimit = 20
    
    #SEEDING
    #VARS
    #print('seeding')
    
    lpRows = model.addVars(rows_sorted[:seed_rows], lb=0, ub=1, vtype=grb.GRB.INTEGER, name='rw')
    lpCols = model.addVars(cols_sorted[:seed_cols], lb=0, ub=1, vtype=grb.GRB.INTEGER, name='cl')
    lpCells = model.addVars([(r, c) for r in rows_sorted[:seed_rows] for c in cols_sorted[:seed_cols]], lb=0, ub=1, vtype=grb.GRB.INTEGER, name='ce')
    #print('size current problem', len(lpRows),len(lpCols))
    #OBJ FCT
    model.setObjective(grb.quicksum([(X_problem[c[0]][c[1]]) * lpCells[c] for c in lpCells]), grb.GRB.MAXIMIZE)
    
    #print()
    #print('Seeding', len(lpRows),len(lpCols))
    #print()
    
    #CONSTRAINTS
    for cell in lpCells:
        model.addConstr(1 - lpRows[cell[0]] >= lpCells[cell], f'{cell}_cr')
        model.addConstr(1 - lpCols[cell[1]] >= lpCells[cell], f'{cell}_cc')
        model.addConstr(1 - lpRows[cell[0]] - lpCols[cell[1]]<= lpCells[cell], f'{cell}_ccr')
        
    model.addConstr(error_rate*grb.quicksum(lpCells) >= grb.quicksum(
            [lpCells[coord]*(1-X_problem[coord[0]][coord[1]]) for coord in lpCells]), 'err_thrshld')
    
    model.optimize()
    
    ##EXTEND BY ROW
    #print('row extend')
    rw = []
    cl = []
    for var in model.getVars():
        if var.X == 0:
            name = var.VarName
            if name[0:2] == 'rw':
                rw += [int(name[3:-1])]
            elif name[0:2] == 'cl':
                cl += [int(name[3:-1])]        
    
    rem_rows = [r for r in rows_sorted if r not in lpRows.keys()]
    rem_rows_sum = X_problem[rem_rows][:,cl].sum(axis=1)
    potential_rows = [r for idx,r in enumerate(rem_rows) if rem_rows_sum[idx]>0.5*len(cl)]
    
    lpRows.update(model.addVars(potential_rows, lb=0, ub=1, vtype=grb.GRB.INTEGER, name='rw'))
    
    new_cells = model.addVars([(r, c) for r in potential_rows for c in cl], lb=0, ub=1, vtype=grb.GRB.INTEGER, name='ce')
    lpCells.update(new_cells)
    
    #print()
    #print('ROW', len(lpRows),len(lpCols))
    #print()
    
    model.setObjective(grb.quicksum([(X_problem[c[0]][c[1]]) * lpCells[c] for c in lpCells]), grb.GRB.MAXIMIZE)
    
    for cell in new_cells:
        model.addConstr(1 - lpRows[cell[0]] >= lpCells[cell], f'{cell}_cr')
        model.addConstr(1 - lpCols[cell[1]] >= lpCells[cell], f'{cell}_cc')
        model.addConstr(1 - lpRows[cell[0]] - lpCols[cell[1]]<= lpCells[cell], f'{cell}_ccr')
        
    model.addConstr(error_rate*grb.quicksum(lpCells) >= grb.quicksum(
            [lpCells[coord]*(1-X_problem[coord[0]][coord[1]]) for coord in lpCells]), 'err_thrshld')
    
    model.optimize()
    
    rw = []
    cl = []
    for var in model.getVars():
        if var.X == 0:
            name = var.VarName
            if name[0:2] == 'rw':
                rw += [int(name[3:-1])]
            elif name[0:2] == 'cl':
                cl += [int(name[3:-1])]
                
    rem_cols = [c for c in cols_sorted if c not in lpCols.keys()]
    rem_cols_sum = X_problem[rw][:,rem_cols].sum(axis=0)
    potential_cols = [c for idx,c in enumerate(rem_cols) if rem_cols_sum[idx]>0.9*len(cl)]
    
    lpCols.update(model.addVars(potential_cols, lb=0, ub=1, vtype=grb.GRB.INTEGER, name='cl'))
    
    new_cells = model.addVars([(r, c) for r in rw for c in potential_cols], lb=0, ub=1, vtype=grb.GRB.INTEGER, name='ce')
    lpCells.update(new_cells)
    
    #print()
    #print('COL', len(lpRows),len(lpCols))
    #print()
    
    model.setObjective(grb.quicksum([(X_problem[c[0]][c[1]]) * lpCells[c] for c in lpCells]), grb.GRB.MAXIMIZE)
    
    for cell in new_cells:
        model.addConstr(1 - lpRows[cell[0]] >= lpCells[cell], f'{cell}_cr')
        model.addConstr(1 - lpCols[cell[1]] >= lpCells[cell], f'{cell}_cc')
        model.addConstr(1 - lpRows[cell[0]] - lpCols[cell[1]]<= lpCells[cell], f'{cell}_ccr')
    
    model.addConstr(error_rate*grb.quicksum(lpCells) >= grb.quicksum(
            [lpCells[coord]*(1-X_problem[coord[0]][coord[1]]) for coord in lpCells]), 'err_thrshld')
    
    
    model.optimize()
    
    rw = []
    cl = []
    for var in model.getVars():
        if var.X == 0:
            name = var.VarName
            if name[0:2] == 'rw':
                rw += [int(name[3:-1])]
            elif name[0:2] == 'cl':
                cl += [int(name[3:-1])]  
    
    # status check
    status = model.Status
    
    if status in (grb.GRB.INF_OR_UNBD, grb.GRB.INFEASIBLE, grb.GRB.UNBOUNDED):
        return (rw,cl,False)
    elif status == grb.GRB.TIME_LIMIT:
        return (rw,cl,True)
    elif status != grb.GRB.OPTIMAL:
        return (rw,cl,False)
               
    return rw,cl,True

def binary_clustering_step(X_matrix, error_rate = 0.025, min_row_quality = 5, min_col_quality = 3):
    
    #Clustering single step
    ###fill the empty space and create an inverse matrix in order to cluster 0
    X_problem = X_matrix
    X_problem_1 = X_problem.copy()
    X_problem_1[X_problem_1 == -1] = 0
    X_problem_0 = X_problem.copy()
    X_problem_0[X_problem_0 == -1] = 1
    X_problem_0 = (X_problem_0-1)*-1
    
    remain_rows = range(X_problem_1.shape[0])
    current_cols = range(X_problem_1.shape[1])
    clustering_1 = True
    status = True
    rw1,rw0 = [],[]
    
    while len(remain_rows)>=min_row_quality and len(current_cols)>=min_col_quality and status:
        if clustering_1:
            rw,cl,status = quasibiclique(X_problem_1[remain_rows][:,current_cols],error_rate)
        else:
            rw,cl,status = quasibiclique(X_problem_0[remain_rows][:,current_cols],error_rate)
             
        rw = [remain_rows[r] for r in rw]
        cl = [current_cols[c] for c in cl]
        
        if len(cl)>0:
            #filter out extremely noisy columns
            if len(rw) == 0 :
                current_cols = []
            else :
                if clustering_1:
                    col_homogeneity = X_problem_1[rw][:,cl].sum(axis=0)/len(rw)
                else:
                    col_homogeneity = X_problem_0[rw][:,cl].sum(axis=0)/len(rw)
                current_cols = [c for idx,c in enumerate(cl) if col_homogeneity[idx] > 5*error_rate]
        else:
            status = False   
            
        if status:
            if clustering_1:
                rw1 = rw1 + [r for r in rw]
            else:
                rw0 = rw0 + [r for r in rw]  
                
        remain_rows = [r for r in remain_rows if r not in rw]
        
        clustering_1 = not clustering_1
    
    return rw1,rw0,current_cols

def biclustering_full_matrix(X_matrix, regions, steps, min_row_quality = 5, min_col_quality = 3,error_rate = 0.025):
    #Iteratively single step cluster through the whole matrix
    print('Clustering', len(regions), 'regions')
    steps_result = steps
    if len(regions)>0:
        print('Clustering region: ', end= "")
    for idx,region in enumerate(regions):   
        remain_cols = region
        
        status = True
        if len(remain_cols)>=min_col_quality:
            while len(remain_cols)>=min_col_quality and status:

                reads1,reads0,cols = binary_clustering_step(X_matrix[:,remain_cols], error_rate = error_rate,min_row_quality = 5, min_col_quality = 3) 
                cols = [remain_cols[c] for c in cols]
                
                if len(cols) == 0:
                    ###not fiding anything of significance so stop to save time 
                    status = False
                else:
                    steps_result.append((reads1,reads0,cols))
                    remain_cols = [c for c in remain_cols if c not in cols]
        print(idx+1, end = " finished ")
    if len(regions)>0:
        print()                           
    return [step for step in steps_result if len(step[0])>0 and len(step[1])>0 and len(step[2])>=min_col_quality]

def post_processing(X_matrix, steps, read_names,distance_thresh = 0.1): 
    #cut the reads into read groups

    ###begin with every reads in the same group
    clusters = [range(len(read_names))]
    
    #go through each steps and seperate them
    for step in steps:
        reads1,reads0,cols = step
        new_clusters = []
        if len(reads0) >0:
            for cluster in clusters:
                clust1 = [c for c in cluster if c in reads1]
                clust0 = [c for c in cluster if c in reads0]
                if len(clust1)>0:
                    new_clusters.append(clust1)
                if len(clust0)>0:
                    new_clusters.append(clust0)
            clusters = new_clusters
    
    ###remove the clusters with less than 5 reads and put them into remaining reads list
    rem_ = []
    big_clusters = []
    for cluster in clusters:
        if len(cluster)<=5:
            rem_ = rem_+ list(cluster)
        else:
            big_clusters.append(cluster)
            
    clusters = big_clusters

    ###calculate the mean vector of each cluster        
    mean_of_clusters = []
    for cluster in clusters:
        mean_of_clusters.append(np.rint(X_matrix[cluster].sum(axis = 0)/len(cluster))) 
    
    ###put remain reads into its closest cluster
    for read in range(len(read_names)):
        is_clustered = False
        for cluster in clusters:
            if read in cluster:
                is_clustered = True
        if not is_clustered:
            rem_.append(read)
    
    if len(rem_) > 0:
        for r in rem_:
            if len(X_matrix[r]) > 0:
                dist = pairwise_distances([X_matrix[r]]+mean_of_clusters, metric = "hamming")[0][1:len(mean_of_clusters)+1]
                if len(dist)>0:
                    idx_most_similar = np.argmin(dist)
                    if dist[idx_most_similar] < 0.1:
                        clusters[idx_most_similar].append(r)
    
    if len(clusters) > 1:
        mean_of_clusters = []
        for cluster in clusters:
            mean_of_clusters.append(np.rint(X_matrix[cluster].sum(axis = 0)/len(cluster)))
        
        agglo_cl = AgglomerativeClustering(n_clusters=None,metric = 'hamming', linkage = 'complete',distance_threshold=distance_thresh)
        agglo_cl.fit(mean_of_clusters)
        labels = agglo_cl.labels_
        new_clusters = {}
        for idx,label in enumerate(labels):
            if label in new_clusters:
                new_clusters[label] = new_clusters[label] + clusters[idx]
            else:
                new_clusters[label] = clusters[idx]    
        
        clusters = new_clusters.values()
    
    ### match indexes to read names
    result_clusters = []
    
    for cluster in clusters:
        result_clusters.append(np.sort(np.array([read_names[r] for r in cluster])))
    
    return result_clusters

def parse_arguments():
    """Parse the input arguments and retrieve the choosen resolution method and
    the instance that must be solve."""
    argparser = ArgumentParser()

    argparser.add_argument(
        '--filepath', dest='filepath', required=True, default='', type=str,
        help='Select the data',
    )

    argparser.add_argument(
        '--error_rate', dest='error_rate', required=False, default=0.025, type=float,
        help='Select the error rate value',
    )

    argparser.add_argument(
        '--out', dest='out', required=True, default='', type=str,
        help='Select the output file',
    )

    argparser.add_argument(
        '--window', dest='window', required=False, default=5000, type=int,
        help='Select the window w',
    )

    arg = argparser.parse_args()


    return (arg.filepath, arg.error_rate, arg.out, arg.window)

if __name__ == '__main__':
    file_path,  error_rate, out, window = parse_arguments()
    start = time.time()
    sol_file = open(out,'w')
    file = ps.AlignmentFile(file_path,'rb')
    contigs = (file.header.to_dict())['SQ']
    if len(contigs) == 0:
        print('ERROR: No contigs found when parsing the BAM file, check the bam file and the indexation of the bam file')
        sys.exit(1)
    for num in range(0,len(contigs)):
        contig_name = contigs[num]['SN']
        contig_length = contigs[num]['LN']

        print(contig_name, contig_length, ' length')
        window = 5000
        filtered_col_threshold = 0.6
        min_row_quality = 5
        min_col_quality = 3
        for start_pos in range(0,contig_length,window):

            if start_pos+window <= contig_length:
                sol_file.write(f'CONTIG\t{contig_name} {start_pos}<->{start_pos+window} \n')

                print(f'Parsing data on contig {contig_name} {start_pos}<->{start_pos+window}')
                dict_of_sus_pos = get_data(file, contig_name,start_pos,start_pos+window)
            else : 
                sol_file.write(f'CONTIG\t{contig_name} {start_pos}<->{contig_length} \n')

                print(f'Parsing data on contig  {contig_name} {start_pos}<->{contig_length}')
                dict_of_sus_pos = get_data(file, contig_name,start_pos,contig_length)

            if len(dict_of_sus_pos) > 0:

                ###create a matrix from the columns
                df = pd.DataFrame(dict_of_sus_pos) 
                # select the reads spanning the whole window
                tmp_idx = df.iloc[:,:len(df.columns)//3].dropna(axis=0, how = 'all')
                df = df.loc[tmp_idx.index,:]
                tmp_idx = df.iloc[:,2*len(df.columns)//3:].dropna(axis=0, how = 'all')
                df = df.loc[tmp_idx.index,:]
                df = df.dropna(axis = 1, thresh = filtered_col_threshold*(len(df.index)))
                reads = list(df.index)

                ###clustering
                X_matrix = df.to_numpy()
                print(X_matrix.shape)
                matrix,regions,steps = pre_processing(X_matrix,min_col_quality)
                steps = biclustering_full_matrix(matrix, regions, steps, min_row_quality, min_col_quality,error_rate=0.025)
                clusters = post_processing(matrix, steps, reads,distance_thresh = 0.05)
                
                reads_ = []
                labels_ = []
                if len(clusters) >1:
                    for idx,cluster in enumerate(clusters):
                        for read in cluster:
                            reads_.append(read)
                            labels_.append(idx)
                    print('Found', len(clusters), 'groups')
                else:
                    print('No haplotypes found')
                    for read in df.index:
                        reads_.append(read)
                        labels_.append(-1)
                        

                for read in reads_:
                    sol_file.write(f'READ\t{read}\n')
                if len(labels_)>0:
                    sol_file.write(f'LABELS\t{labels_[0]}')
                    for i in range(1,len(labels_)):
                        sol_file.write(f',{labels_[i]}')
                
                sol_file.write(f'\n')   
                end = time.time()
                print('Elapsed time', end - start)
            else:
                print('No suspicious positions found')


    sol_file.close()    