#ifndef GRAPHUNZIP_H
#define GRAPHUNZIP_H

#include <iostream>
#include <vector>
#include <string>
#include <mutex>
#include <set>
#include "robin_hood.h"

namespace std {
    template <>
    class hash<std::pair<int, int>> {
    public:
        size_t operator()(const std::pair<int, bool>& pair) const {
            return std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second);
        }
    };
}

#endif