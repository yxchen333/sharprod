#include "random.h"

std::mt19937_64 gen(std::random_device{}());
std::uniform_real_distribution<double> dis(0.0, 1.0);
