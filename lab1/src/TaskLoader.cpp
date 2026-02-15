#include "../include/TaskLoader.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace
{

	struct VarBoundInfo
	{
		int var_index;
		BoundType type;
	};

	VarBoundInfo parseVarBound(const string &line, int n_vars)
	{
		istringstream iss(line);
		string var_part, op, val;
		iss >> var_part >> op >> val;

		if (var_part.empty() || var_part[0] != 'x')
			throw runtime_error("Invalid variable format: " + var_part);
		int var_num = stoi(var_part.substr(1));

		if (op != ">=" || val != "0")
			throw runtime_error("Unsupported bound: " + op + " " + val);

		if (var_num < 1 || var_num > n_vars)
			throw runtime_error("Variable index out of range: x" + to_string(var_num));

		return {var_num, BoundType::NOT_NEGATIVE};
	}

} // namespace

vector<string> LPProblemParser::split(const string &s)
{
	vector<string> tokens;
	istringstream iss(s);
	string token;
	while (iss >> token)
		tokens.push_back(token);
	return tokens;
}

InequalityType LPProblemParser::parseInequalityType(const string &op)
{
	if (op == "==")
		return InequalityType::EQUAL;
	if (op == "<=")
		return InequalityType::LESS_EQUAL;
	if (op == ">=")
		return InequalityType::GREATER_EQUAL;
	throw invalid_argument("Invalid operator: " + op);
}
int isSubstring(string s1, string s2)
{
	int M = s1.length();
	int N = s2.length();

	/* A loop to slide pat[] one
	   by one */
	for (int i = 0; i <= N - M; i++)
	{
		int j;

		/* For current index i, check for
		   pattern match */
		for (j = 0; j < M; j++)
			if (s2[i + j] != s1[j])
				break;

		if (j == M)
			return i;
	}
	return -1;
}

LPProblemSlack *LPProblemParser::parse(const string &filename)
{
	ifstream file(filename);
	if (!file)
		throw runtime_error("Cannot open file: " + filename);

	string line;

	// Read n and m
	getline(file, line);
	auto tokens = split(line);
	if (tokens.size() != 2)
		throw runtime_error("Invalid n m format");
	int n = stoi(tokens[0]), m = stoi(tokens[1]);

	// Read objective type
	getline(file, line);
	ObjectiveType objType = (isSubstring("min",line) != -1) ? ObjectiveType::MINIMIZE : ObjectiveType::MAXIMIZE;

	// Read objective coefficients
	getline(file, line);
	tokens = split(line);
	if (tokens.size() != n)
		throw runtime_error("Invalid objective coefficients count");
	vector<double> c(n);
	transform(tokens.begin(), tokens.end(), c.begin(), [](const string &s)
			  { return stod(s); });

	// Read constraints
	vector<Constraint> constraints;
	int eqCount = 0, leCount = 0, geCount = 0;
	for (int i = 0; i < m; ++i)
	{
		getline(file, line);
		tokens = split(line);
		if (tokens.size() != n + 2)
			throw runtime_error("Invalid constraint format");

		vector<double> coeffs(n);
		transform(tokens.begin(), tokens.begin() + n, coeffs.begin(), [](const string &s)
				  { return stod(s); });

		InequalityType type = parseInequalityType(tokens[n]);
		double rhs = stod(tokens[n + 1]);

		// Создаём объект основного Constraint
		Constraint c;
		c.coefficients = move(coeffs);
		c.b = rhs; // Используем поле b вместо rhs
		c.type = type;

		constraints.push_back(c);

		switch (type)
		{
		case InequalityType::EQUAL:
			eqCount++;
			break;
		case InequalityType::LESS_EQUAL:
			leCount++;
			break;
		case InequalityType::GREATER_EQUAL:
			geCount++;
			break;
		}
	}


	// Read variable bounds
	vector<VarBoundInfo> varBounds;
	while (getline(file, line))
	{
		if (line.empty())
			continue;
		varBounds.push_back(parseVarBound(line, n));
	}


	// Create problem
	LPProblemSlack *problem = new LPProblemSlack(n);
	problem->set_objective(c, objType);
	for (const auto &con : constraints)
	{
		problem->add_constraint(con); // Теперь типы совпадают
	}

	for (const auto &vb : varBounds)
		problem->add_var_bound({vb.var_index, vb.type});

	return problem;
}