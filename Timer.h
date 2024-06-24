//
// Created by Stepan Didurenko on 24.06.2024.
//

#ifndef CPP_LAB_15_16_TIMER_H
#define CPP_LAB_15_16_TIMER_H

#include <chrono>

class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
public:
    Timer() = default;

    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    std::chrono::milliseconds get() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_time);
    }
};

#endif //CPP_LAB_15_16_TIMER_H
