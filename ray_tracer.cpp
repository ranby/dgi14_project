#include <iostream>
//#include <glm/glm.hpp>
#include </Users/galgazur/Downloads/CgLab1/glm/glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
using glm::vec3;
using glm::mat3;

// ----------------------------------------------------------------------------
// STRUCTS

struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
};

struct Lense
{
	vec3 center;
	float radius;
	vec3 normal;
	float focalLength;
	float refractiveIndex;
};

struct LenseIntersection
{
	vec3 position;
	int lenseIndex;
	float distance;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 300;
const int SCREEN_HEIGHT = 300;
SDL_Surface* screen;
int t;
float PI = 3.14159f;

float focalLength = SCREEN_HEIGHT / 2;
vec3 cameraPos(0, 0, -2.8);
vector<Triangle> triangles;
mat3 R;
float yaw;

vec3 lightPos(0, -0.5, -0.7);
vec3 lightColor = 14.f * vec3(1, 1, 1);
vec3 indirectLight = 0.5f*vec3(1, 1, 1);

vector<Lense> lenses;
float defaultRefractiveIndex = 1.f;

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void CalculateRotation();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
bool GetIntersectedTriangleColor(vec3 startPosition, vec3 dir, vec3& intersectionColor, Intersection& isn);
bool Intersects(vec3 x);
vec3 DirectLight(const Intersection& i);

void SetupLenses();
bool IntersectsLense(vec3 start, vec3 dir, LenseIntersection& intersection, int previousIndex);
void calculateRefraction(vec3 dirIn, vec3 lensePointIn, int lenseIndex, vec3& lensePointOut, vec3& dirOut);
void calculateRefractionVector(Lense lense, vec3 dirIn, vec3 pointIn, float mediumIn, float mediumOut, vec3& dirOut);
void calculateReflection(vec3 dirIn, vec3 lensePoint, Lense lense, vec3& dirOut);
void findPerpendicular(vec3 aVector, vec3& perpendicularVector);
vec3 rotateVector(vec3 vector, float pitch, float yaw);


// ----------------------------------------------------------------------------
// CODE

void noop(int x, int y){
	if (x == 35){
		if (y == 50){
			int i = 0;
		}
	}
}

void findPerpendicular(vec3 aVector, vec3& perpendicularVector){
	srand(time(NULL));

	vec3 bVector;
	bVector = vec3(rand() % 50, rand() % 50, rand() % 50);
	bVector = glm::normalize(bVector);

	if (glm::dot(aVector, bVector) == 0){
		bVector.x = bVector.x + 10;
	}

	perpendicularVector = glm::cross(aVector, bVector);
}

int main(int argc, char* argv[])
{
	yaw = 0;
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);

	SetupLenses();

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	vec3 right(R[0][0], R[0][1], R[0][2]);
	vec3 down(R[1][0], R[1][1], R[1][2]);
	vec3 forward(R[2][0], R[2][1], R[2][2]);

	Uint8* keystate = SDL_GetKeyState(0);
	if (keystate[SDLK_UP])
	{
		cameraPos = cameraPos + forward / float(4);
	}
	if (keystate[SDLK_DOWN])
	{
		cameraPos = cameraPos - forward / float(4);
	}
	if (keystate[SDLK_LEFT])
	{
		yaw = yaw - float(PI / 8);
		CalculateRotation();
	}
	if (keystate[SDLK_RIGHT])
	{
		yaw = yaw + float(PI / 8);
		CalculateRotation();
	}
	if (keystate[SDLK_w]) {
		lightPos += forward / float(4);
	}
	if (keystate[SDLK_s]) {
		lightPos -= forward / float(4);
	}
	if (keystate[SDLK_a]) {
		lightPos -= right / float(4);
	}
	if (keystate[SDLK_d]) {
		lightPos += right / float(4);
	}
	if (keystate[SDLK_q]) {
		lightPos -= down / float(4);
	}
	if (keystate[SDLK_e]) {
		lightPos += down / float(4);
	}

}

void CalculateRotation() {
	vec3 col1(cos(yaw), 0, sin(yaw));
	vec3 col2(0, 1, 0);
	vec3 col3(-sin(yaw), 0, cos(yaw));
	R = mat3(col1, col2, col3);
}

void Draw()
{
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	Intersection isn1;
	Intersection isn2;
	isn2.distance = INT_MAX;
	for (int y = 0; y<SCREEN_HEIGHT; ++y)
	{
		for (int x = 0; x<SCREEN_WIDTH; ++x)
		{
			noop(x, y);

			vec3 triangleColor(0, 0, 0);
			vec3 d(x - (SCREEN_WIDTH / 2), y - (SCREEN_HEIGHT / 2), focalLength);
			d = d*R;

			LenseIntersection lenseIntersection;
			lenseIntersection.position.z = INT_MAX;
			vec3 lenseColor = vec3(0, 0, 0);
			if (IntersectsLense(cameraPos, d, lenseIntersection, -1)) {

				vec3 pointOut;
				vec3 dirOut;
				calculateRefraction(d, lenseIntersection.position, lenseIntersection.lenseIndex, pointOut, dirOut);

				GetIntersectedTriangleColor(pointOut, dirOut, lenseColor, isn1);
			}
			float lenseDistance = glm::length(lenseIntersection.position - cameraPos);

			GetIntersectedTriangleColor(cameraPos, d, triangleColor, isn2);

			vec3 color;
			if (lenseDistance < isn2.distance) {
				color = lenseColor;
			}
			else {
				color = triangleColor;
			}
			PutPixelSDL(screen, x, y, color);
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

bool GetIntersectedTriangleColor(vec3 startPosition, vec3 dir, vec3& intersectionColor, Intersection& isn) {
	if (ClosestIntersection(startPosition, dir, triangles, isn)){
		vec3 illumination = DirectLight(isn);
		intersectionColor = triangles[isn.triangleIndex].color*(illumination + indirectLight);
		return true;
	}
	return false;
}


bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection)
{
	Intersection tempI;

	float m = std::numeric_limits<float>::max();
	float shortestDistance = m;

	for (size_t i = 0; i<triangles.size(); i++)
	{
		Triangle triangle = triangles[i];
		vec3 v0 = triangle.v0;
		vec3 v1 = triangle.v1;
		vec3 v2 = triangle.v2;
		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = start - v0;
		mat3 A(-dir, e1, e2);
		vec3 x = glm::inverse(A) * b;


		tempI.position = v0 + x.y*e1 + x.z*e2;
		tempI.distance = glm::length(tempI.position - start);
		tempI.triangleIndex = i;

		if (Intersects(x) && tempI.distance > 0.0001f)
		{
			if (shortestDistance == m || tempI.distance <= shortestDistance)
			{
				shortestDistance = tempI.distance;
				closestIntersection = tempI;
			}
		}

	}
	if (shortestDistance != m) {
		return true;
	}
	else {
		return false;
	}
}


bool Intersects(vec3 x)
{
	if (x.y<0)
		return false;

	if (x.z<0)
		return false;

	if (x.y + x.z>1)
		return false;

	if (x.x<0)
		return false;

	return true;
}


vec3 DirectLight(const Intersection& i){

	vec3 rVector = lightPos - i.position;
	Triangle triangle = triangles[i.triangleIndex];
	vec3 normal = triangle.normal;

	Intersection isn;
	if (ClosestIntersection(i.position, rVector, triangles, isn))
	{
		if ((glm::abs(isn.distance) <= glm::abs(glm::length(rVector))))
		{
			return vec3(0, 0, 0);
		}
	}

	float rLength = glm::length(rVector);
	float dotProd = glm::dot(rVector / rLength, normal);
	vec3 illumination = (lightColor*fmax(dotProd, 0.0f)) / (4 * rLength*rLength*(PI));

	return illumination;
}

void SetupLenses() {
	lenses = vector<Lense>(2);

	Lense lense;
	lense.center = vec3(0, 0, -1.3f);
	lense.radius = 0.7f;
	lense.normal = glm::normalize(vec3(0, 0, -1));
	lense.focalLength = -3.5f;
	lense.refractiveIndex = 1.5f;
	lenses[0] = lense;

	lense.center = vec3(0, 0, -1.3f);
	lense.radius = 0.7f;
	lense.normal = glm::normalize(vec3(0, 0, 1));
	lense.focalLength = -3.f;
	lense.refractiveIndex = 1.5f;
	lenses[1] = lense;
}

bool IntersectsLense(vec3 start, vec3 dir, LenseIntersection& intersection, int previousIndex) {
	if (lenses.size() <= 0) {
		return false;
	}
	dir = glm::normalize(dir);
	Lense lense;
	for (int i = 0; i < lenses.size(); ++i) {
		if (i != previousIndex) {
			lense = lenses[i];

			vec3 focalpoint = (-lense.focalLength)*lense.normal + lense.center;

			vec3 ray = start - focalpoint;
			float dotProd = glm::dot(dir, ray);
			float length = glm::length(ray);

			vec3 np;
			findPerpendicular(lense.normal, np);
			np = glm::normalize(np);
			np = lense.center + (np*lense.radius);

			float sphereRadius = glm::length(focalpoint - np);

			float square = (dotProd * dotProd) - (length * length) + (sphereRadius * sphereRadius);

			float d1, d2;

			if (square == 0){
				d1 = -dotProd;
				d2 = INT_MAX;
			}
			else if (square > 0){
				d1 = -dotProd + sqrt(square);
				d2 = -dotProd - sqrt(square);
			}
			else {
				d1 = INT_MAX;
				d2 = INT_MAX;

			}

			vec3 intersection1 = start + d1*dir;
			vec3 intersection2 = start + d2*dir;

			if (d1 != INT_MAX){
				if ((glm::dot(intersection1 - lense.center, lense.normal) > 0) && (glm::length(intersection1 - start) > 0)){
					intersection.position = intersection1;
					intersection.lenseIndex = i;
					return true;
				}
				else if (((d2 != INT_MAX) && glm::dot((intersection2 - lense.center), lense.normal) > 0) && (glm::length(intersection2 - start) > 0)){
					intersection.position = intersection2;
					intersection.lenseIndex = i;
					return true;
				}
			}
		}
	}
	return false;
}

void calculateRefraction(vec3 dirIn, vec3 lensePointIn, int lenseIndex, vec3& lensePointOut, vec3& dirOut) {
	vec3 insideRefractionDir;
	Lense lense = lenses[lenseIndex];
	calculateRefractionVector(lense, dirIn, lensePointIn, defaultRefractiveIndex, lense.refractiveIndex, insideRefractionDir);

	LenseIntersection outIntersection;
	if (IntersectsLense(lensePointIn, insideRefractionDir, outIntersection, lenseIndex)) {
		Lense outLense = lenses[outIntersection.lenseIndex];
		calculateRefractionVector(outLense, insideRefractionDir, outIntersection.position, outLense.refractiveIndex, defaultRefractiveIndex, dirOut);

		lensePointOut = outIntersection.position;
	}
	else {
		//much error, such bad, wow
	}
}

void calculateRefractionVector(Lense lense, vec3 dirIn, vec3 pointIn, float mediumIn, float mediumOut, vec3& dirOut) {
	//finding normal for the point on the lense where the ray hits
	vec3 focalPoint = (-lense.focalLength)*lense.normal + lense.center;
	vec3 normalVector = glm::normalize(pointIn - focalPoint);

	float flipDot = glm::dot(normalVector, dirIn); // check if ray is on its way in or out of the lense
	if (flipDot < 0){
		normalVector = glm::normalize(focalPoint - pointIn);
	}

	//calculating the angle between incoming ray and previous found normal

	if (dirIn.x == 0)
		dirIn.x += 0.00000001f;
	if (dirIn.y == 0)
		dirIn.y += 0.00000001f;
	if (dirIn.z == 0)
		dirIn.z += 0.00000001f;
	if (normalVector.x == 0)
		normalVector.x += 0.00000001f;
	if (normalVector.y == 0)
		normalVector.y += 0.00000001f;
	if (normalVector.z == 0)
		normalVector.z += 0.00000001f;
	
	//float pitch = glm::radians(90.f) - glm::atan(dirIn.x / dirIn.y) - glm::atan(normalVector.y / normalVector.x);
	//float yaw = glm::radians(90.f) - glm::atan(dirIn.x / dirIn.z) - glm::atan(normalVector.z / normalVector.x);
	float pitch1 = glm::atan(dirIn.y / dirIn.z) - glm::atan(normalVector.y / normalVector.z);
	float yaw1 = glm::atan(dirIn.x / dirIn.z) - glm::atan(normalVector.x / normalVector.z);

	//Snell's law
	float pitch2 = glm::asin(mediumIn * glm::sin(pitch1) / mediumOut);
	float yaw2 = glm::asin(mediumIn * glm::sin(yaw1) / mediumOut);

	vec3 n2;
	n2 = normalVector;

//	vec3 pitchedVector = vec3(dirIn.x, (dirIn.y * glm::cos(pitch)) + (dirIn.x * glm::sin(pitch)), (dirIn.z * glm::cos(pitch)) - (dirIn.x * glm::sin(pitch)));
//	vec3 yawedVector = vec3((pitchedVector.x * glm::cos(yaw)) - (pitchedVector.z * glm::sin(yaw)), pitchedVector.y, (pitchedVector.z * glm::cos(yaw)) - (pitchedVector.x * glm::sin(yaw)));

	vec3 rotatedVector = rotateVector(dirIn, pitch2-pitch1, yaw1-yaw2);

	dirOut = rotatedVector;
}

vec3 rotateVector(vec3 vector, float pitch, float yaw) {
	mat3 matPitch = mat3(vec3(1, 0, 0), vec3(0, cos(pitch), -sin(pitch)), vec3(0, sin(pitch), cos(pitch)));
	mat3 matYaw = mat3(vec3(cos(yaw), 0, sin(yaw)), vec3(0, 1, 0), vec3(-sin(yaw), 0, cos(yaw)));
	mat3 mat = matPitch * matYaw;
	vec3 rotated = mat * vector;
	return rotated;
}

void calculateReflection(vec3 dirIn, vec3 lensePoint, vec3& dirOut) {
	//TODO
}
