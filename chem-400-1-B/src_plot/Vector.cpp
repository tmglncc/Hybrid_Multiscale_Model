#include "Vector.hpp"

Vector3::Vector3(double x, double y, double z): x(x), y(y), z(z) {}

Vector3 Vector3::operator+(const Vector3& v) {
	return Vector3(this->x + v.x, this->y + v.y, this->z + v.z);
}

Vector3 Vector3::operator-(const Vector3& v) {
	return Vector3(this->x - v.x, this->y - v.y, this->z - v.z);
}

Vector3 Vector3::operator+(const double c) {
	return Vector3(this->x + c, this->y + c, this->z + c);
}

Vector3 Vector3::operator-(const double c) {
	return Vector3(this->x - c, this->y - c, this->z - c);
}

void Vector3::operator+=(const Vector3& v) {
	this->x += v.x;
	this->y += v.y;
	this->z += v.z;
}

void Vector3::operator-=(const Vector3& v) {
	this->x -= v.x;
	this->y -= v.y;
	this->z -= v.z;
}

void Vector3::operator+=(const double c) {
	this->x += c;
	this->y += c;
	this->z += c;
}

void Vector3::operator-=(const double c) {
	this->x -= c;
	this->y -= c;
	this->z -= c;
}

Vector3 Vector3::operator*(const Vector3& v) {
	return Vector3(this->x * v.x, this->y * v.y, this->z * v.z);
}

Vector3 Vector3::operator/(const Vector3& v) {
	return Vector3(this->x / v.x, this->y / v.y, this->z / v.z);
}

Vector3 Vector3::operator*(const double c) {
	return Vector3(this->x * c, this->y * c, this->z * c);
}

Vector3 Vector3::operator/(const double c) {
	return Vector3(this->x / c, this->y / c, this->z / c);
}

void Vector3::operator*=(const Vector3& v) {
	this->x *= v.x;
	this->y *= v.y;
	this->z *= v.z;
}

void Vector3::operator/=(const Vector3& v) {
	this->x /= v.x;
	this->y /= v.y;
	this->z /= v.z;
}

void Vector3::operator*=(const double c) {
	this->x *= c;
	this->y *= c;
	this->z *= c;
}

void Vector3::operator/=(const double c) {
	this->x /= c;
	this->y /= c;
	this->z /= c;
}

bool Vector3::operator==(const Vector3& vec) {
	return this->x == vec.x && this->y == vec.y && this->z == vec.z;
}

bool Vector3::operator!=(const Vector3& vec) {
	return !(this->x == vec.x && this->y == vec.y && this->z == vec.z);
}

std::string Vector3::to_string() {
	return "(" + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + ")";
}

Vector3i::Vector3i(int x, int y, int z): x(x), y(y), z(z) {}

Vector3i::Vector3i(Vector3 vec) {
	this->x = (int) vec.x;
	this->y = (int) vec.y;
	this->z = (int) vec.z;
}

Vector3i Vector3i::operator+(const Vector3i& v) {
	return Vector3i(this->x + v.x, this->y + v.y, this->z + v.z);
}

Vector3i Vector3i::operator-(const Vector3i& v) {
	return Vector3i(this->x - v.x, this->y - v.y, this->z - v.z);
}

Vector3i Vector3i::operator+(const int c) {
	return Vector3i(this->x + c, this->y + c, this->z + c);
}

Vector3i Vector3i::operator-(const int c) {
	return Vector3i(this->x - c, this->y - c, this->z - c);
}

void Vector3i::operator+=(const Vector3i& v) {
	this->x += v.x;
	this->y += v.y;
	this->z += v.z;
}

void Vector3i::operator-=(const Vector3i& v) {
	this->x -= v.x;
	this->y -= v.y;
	this->z -= v.z;
}

void Vector3i::operator+=(const int c) {
	this->x += c;
	this->y += c;
	this->z += c;
}

void Vector3i::operator-=(const int c) {
	this->x -= c;
	this->y -= c;
	this->z -= c;
}

Vector3i Vector3i::operator*(const Vector3i& v) {
	return Vector3i(this->x * v.x, this->y * v.y, this->z * v.z);
}

Vector3i Vector3i::operator/(const Vector3i& v) {
	return Vector3i(this->x / v.x, this->y / v.y, this->z / v.z);
}

Vector3i Vector3i::operator*(const int c) {
	return Vector3i(this->x * c, this->y * c, this->z * c);
}

Vector3i Vector3i::operator/(const int c) {
	return Vector3i(this->x / c, this->y / c, this->z / c);
}

void Vector3i::operator*=(const Vector3i& v) {
	this->x *= v.x;
	this->y *= v.y;
	this->z *= v.z;
}

void Vector3i::operator/=(const Vector3i& v) {
	this->x /= v.x;
	this->y /= v.y;
	this->z /= v.z;
}

void Vector3i::operator*=(const int c) {
	this->x *= c;
	this->y *= c;
	this->z *= c;
}

void Vector3i::operator/=(const int c) {
	this->x /= c;
	this->y /= c;
	this->z /= c;
}

bool Vector3i::operator==(const Vector3i& vec) {
	return this->x == vec.x && this->y == vec.y && this->z == vec.z;
}

bool Vector3i::operator!=(const Vector3i& vec) {
	return !(this->x == vec.x && this->y == vec.y && this->z == vec.z);
}

std::string Vector3i::to_string() {
	return "(" + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + ")";
}