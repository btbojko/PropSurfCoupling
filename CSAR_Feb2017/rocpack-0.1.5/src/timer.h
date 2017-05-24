#ifndef PACK_TIMER_H
#define PACK_TIMER_H

#include <chrono>

using namespace std::chrono;

template <class unit = nanoseconds>
class Timer {
private:
	using clock = high_resolution_clock;
	clock::time_point t_start;

public:
	Timer() : t_start() { start(); }

	void reset() { start(); }

	void start() { t_start = clock::now(); }

	double elapsed()
	{
		return duration_cast<unit>(clock::now() - t_start).count();
	}
};

class cycles {};

template <>
class Timer<cycles> {
private:
	unsigned long long int t_start;

public:
	Timer() { start(); }

	void reset() { start(); }

	void start() { t_start = get_cycle_count(); }

	double Elapsed()
	{
		return static_cast<double>(get_cycle_count() - t_start);
	}

private:
	inline __attribute__((always_inline))
	unsigned long long int get_cycle_count()
	{
		unsigned int hi, lo;
		asm volatile("cpuid\n\t" "rdtsc" : "=a"(lo), "=d"(hi));
		return ((unsigned long long)lo)|(((unsigned long long)hi) << 32);
	}
};

#endif
