#pragma once

#include "SimplexEngine.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

class LPProblemParser
{
public:
	static std::unique_ptr<LPProblem> parse(const std::string &filename);

private:
	static std::vector<std::string> split(const std::string &s);
	static InequalityType parseInequalityType(const std::string &op);
};