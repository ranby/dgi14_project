#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

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

const int SCREEN_WIDTH = 100;
const int SCREEN_HEIGHT = 100;
SDL_Surface* screen;
int t;
float PI = 3.14159f;

float focalLength = SCREEN_HEIGHT / 2;
vec3 cameraPos(0, 0, -2);
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
bool Intersects(vec3 x);
vec3 DirectLight(const Intersection& i);

void SetupLenses();
bool IntersectsLense(vec3 start, vec3 dir, LenseIntersection& intersection);
void calculateRefraction(vec3 dirIn, vec3 lensePointIn, Lense lense, vec3& lensePointOut, vec3& dirOut);
void calculateReflection(vec3 dirIn, vec3 lensePoint, Lense lense, vec3& dirOut);



// ----------------------------------------------------------------------------
// CODE

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

	Intersection isn;
	for (int y = 0; y<SCREEN_HEIGHT; ++y)
	{
		for (int x = 0; x<SCREEN_WIDTH; ++x)
		{
			vec3 triangleColor(0, 0, 0);
			vec3 d(x - (SCREEN_WIDTH / 2), y - (SCREEN_HEIGHT / 2), focalLength);
			d = d*R;

			LenseIntersection intersection;
			if (IntersectsLense(cameraPos, d, intersection)) {



				//TODO
			}
			float lenseDistance = glm::length(intersection.position - cameraPos);
			vec3 lenseColor = vec3(0, 0, 0); //TODO

			if (ClosestIntersection(cameraPos, d, triangles, isn)){
				vec3 illumination = DirectLight(isn);
				triangleColor = triangles[isn.triangleIndex].color*(illumination + indirectLight);
			}
			
			vec3 color;
			if (lenseDistance < isn.distance) {
				color = lenseColor;
			} else {
				color = triangleColor;
			}
			PutPixelSDL(screen, x, y, color);
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
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

//DEPRECATED MOTHAFUCKA!!
void SetupLenses() {
	Lense lense;
	lense.center = vec3(0, 0, 0);
	lense.radius = 2.f;
	lense.normal = vec3(0, 0, 1);
	lense.focalLength = 3.5f;
	//lense.backFocalLength = 3.5f;
	lense.refractiveIndex = 1.5f;
	lenses[0] = lense;
}

bool IntersectsLense(vec3 start, vec3 dir, LenseIntersection& intersection) {
	Lense lense = lenses[0];

	float dotProd = glm::dot(dir, (start - lense.center));
	float length = glm::length(start - lense.center);
	float square = (dotProd * dotProd) - (length * length) + (lense.radius * lense.radius);

	float d1, d2;

	if (square == 0){
		d1 = -glm::dot(dir, (start - lense.center));
		d2 = INT_MAX;

	}
	else if (square>0){
		d1 = -glm::dot(dir, (start - lense.radius)) + sqrt(square);
		d2 = -glm::dot(dir, (start - lense.radius)) - sqrt(square);
	}
	else {
		d1 = INT_MAX;
		d2 = INT_MAX;

	}

	vec3 intersection1 = start + d1*dir;
	vec3 intersection2 = start + d2*dir;


	if (d1 == INT_MAX){
		return false;
	}
	else{
		if (glm::dot(intersection1 - lense.center, lense.normal) > 0){
			intersection.position = intersection1;
			return true;
		}
		else if ((d2 != INT_MAX) && glm::dot((intersection2 - lense.center), lense.normal) > 0){
			intersection.position = intersection2;
			return true;
		}
	}
}

void calculateRefraction(vec3 dirIn, vec3 lensePointIn, Lense lense, vec3& lensePointOut, vec3& dirOut) {
	vec3 insideRefractionDir;
	calculateRefractionVector(lense, dirIn, lensePointIn, defaultRefractiveIndex, lense.refractiveIndex, insideRefractionDir);

	LenseIntersection outIntersection;
	if (IntersectsLense(lensePointIn, insideRefractionDir, outIntersection)) {
		calculateRefractionVector(lense, insideRefractionDir, outIntersection.position, lense.refractiveIndex, defaultRefractiveIndex, dirOut);

		lensePointOut = outIntersection.position;
	}
	else {
		//much error, such bad, wow
	}
}

void calculateRefractionVector(Lense lense, vec3 dirIn, vec3 pointIn, float mediumIn, float mediumOut, vec3& dirOut) {
	//finding normal for the point on the lense where the ray hits
	vec3 focalPoint = (lense.normal * lense.focalLength) + lense.center;
	vec3 normalVector = glm::normalize(pointIn - focalPoint);

	//calculating the angle between incoming ray and previous found normal
	float pitch = glm::radians(90) - glm::atan(dirIn.x / dirIn.y) - glm::atan(normalVector.y / normalVector.x);
	float yaw = glm::radians(90) - glm::atan(dirIn.x / dirIn.z) - glm::atan(normalVector.z / normalVector.x);

	//Snell's law
	pitch = glm::asin(mediumIn * glm::sin(pitch) / mediumOut);
	yaw = glm::asin(mediumIn * glm::sin(yaw) / mediumOut);

	vec3 n2 = -normalVector;
	vec3 pitchedVector = vec3(n2.x * glm::cos(pitch) - n2.y * glm::sin(pitch), n2.y * glm::cos(pitch) + n2.x * glm::sin(pitch), n2.z);
	vec3 yawedVector = vec3(pitchedVector.x * glm::cos(yaw) - pitchedVector.z * glm::cos(yaw), pitchedVector.y, pitchedVector.z * glm::cos(yaw) - pitchedVector.x * glm::sin(yaw));

	dirOut = yawedVector;
}

void calculateReflection(vec3 dirIn, vec3 lensePoint, vec3& dirOut) {
	//TODO
}
