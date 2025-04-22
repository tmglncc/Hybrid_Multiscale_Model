// #pragma once
#ifndef VECTOR
#define VECTOR

#include <string>

class Vector3 {
public:
	double x, y, z;

	Vector3(double x = 0.0, double y = 0.0, double z = 0.0);
	Vector3 operator+(const Vector3& v);
	Vector3 operator-(const Vector3& v);
	Vector3 operator+(const double c);
	Vector3 operator-(const double c);
	void operator+=(const Vector3& v);
	void operator-=(const Vector3& v);
	void operator+=(const double c);
	void operator-=(const double c);
	Vector3 operator*(const Vector3& v);
	Vector3 operator/(const Vector3& v);
	Vector3 operator*(const double c);
	Vector3 operator/(const double c);
	void operator*=(const Vector3& v);
	void operator/=(const Vector3& v);
	void operator*=(const double c);
	void operator/=(const double c);
	bool operator==(const Vector3& vec);
	bool operator!=(const Vector3& vec);
	std::string to_string();
};

class Vector3i {
public:
	int x, y, z;

	Vector3i(int x = 0.0f, int y = 0.0f, int z = 0.0f);
	Vector3i(Vector3 vec);
	Vector3i operator+(const Vector3i& v);
	Vector3i operator-(const Vector3i& v);
	Vector3i operator+(const int c);
	Vector3i operator-(const int c);
	void operator+=(const Vector3i& v);
	void operator-=(const Vector3i& v);
	void operator+=(const int c);
	void operator-=(const int c);
	Vector3i operator*(const Vector3i& v);
	Vector3i operator/(const Vector3i& v);
	Vector3i operator*(const int c);
	Vector3i operator/(const int c);
	void operator*=(const Vector3i& v);
	void operator/=(const Vector3i& v);
	void operator*=(const int c);
	void operator/=(const int c);
	bool operator==(const Vector3i& vec);
	bool operator!=(const Vector3i& vec);
	std::string to_string();
};

#endif
