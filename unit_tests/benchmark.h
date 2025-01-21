#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>

#define BENCHMARK(function) { auto start = std::chrono::high_resolution_clock::now();\
	function;\
	auto end = std::chrono::high_resolution_clock::now();\
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();\
	std::cout << ">..benchmark@ " << std::left << std::setw(30) << #function << "| excution time : " << duration << " ms :D" << std::endl; }



#define BENCHMARK_BEGIN(token) auto start_##token = std::chrono::high_resolution_clock::now()
#define BENCHMARK_END(token)auto end_##token = std::chrono::high_resolution_clock::now();\
	auto duration_##token = std::chrono::duration_cast<std::chrono::milliseconds>(end_##token - start_##token).count();\
	std::cout << ">..benchmark@ " << std::left << std::setw(30) << #token << "| excution time : " << duration_##token << " ms :D" << std::endl