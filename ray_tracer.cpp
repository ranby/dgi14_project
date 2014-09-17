#include <iostream>
#include <string>
#include <sstream>

#include <glm/glm.hpp>
//#include </Users/galgazur/Downloads/CgLab1/glm/glm/glm.hpp>
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
	int index;
};

struct LenseIntersection
{
	vec3 position;
	int lenseIndex;
	float distance;
};

struct LenseNoiseMap
{
	float maxangle;
	vector<vector<vec3>> noisematrix;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 400;
const int SCREEN_HEIGHT = 400;
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

float focals1[] = { 0.001, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };


float focals2[] = { 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0 };

vector<Lense> lenses;
vector<LenseNoiseMap> lensenoises;

vec3 lenseCenterStart = vec3(0, 0, -1.5f);
float defaultRefractiveIndex = 1.f;

float maxangle;
int noisecells = 150;

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Update2(int index);
void CalculateRotation();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
bool GetIntersectedTriangleColor(vec3 startPosition, vec3 dir, vec3& intersectionColor, Intersection& isn);
bool Intersects(vec3 x);
vec3 DirectLight(const Intersection& i);

void SetupLenses(vec3 center);
bool IntersectsLense(vec3 start, vec3 dir, LenseIntersection& intersection, int previousIndex);
void calculateRefraction(vec3 dirIn, vec3 lensePointIn, int lenseIndex, vec3& lensePointOut, vec3& dirOut);
void calculateRefractionVector(Lense lense, vec3 dirIn, vec3 pointIn, float mediumIn, float mediumOut, vec3& dirOut);
void calculateReflection(vec3 dirIn, vec3 lensePoint, Lense lense, vec3& dirOut);
void findPerpendicular(vec3 aVector, vec3& perpendicularVector);
vec3 rotateVector(vec3 vector, float pitch, float yaw);
void SetupLenseNoises(vector<Lense> lenses);
vec3 calculateNoise(int lenseindex, float anglex, float angley);
float calculateMaxAngle(Lense lense);
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

	SetupLenses(lenseCenterStart);
	//for (int xcv = 0; xcv<200; xcv++) {
	for (int pic = 199; pic<200; pic++) {
		SetupLenseNoises(lenses);
		Update2(pic);
		Draw();

		const char * filename;// = "'~/Desktop/lense_screenshot";

		std::stringstream sstm;
		sstm << "C:\\Users\\Erik\\Documents\\Visual Studio 2013\\Projects\\dgi_project\\dgi_project\\screen" << pic << "-" << ".bmp";
		std:string sstmstr = sstm.str();
		filename = sstmstr.c_str();

		cout << SDL_SaveBMP(screen, filename) << ":" << pic << "%  " << filename << "\n";
		Update2(pic);
	}
	//}
	//	while (NoQuitMessageSDL())
	//	{
	//		Update();
	//		Draw();
	//	}

	//SDL_SaveBMP(screen, "screenshot.bmp");

	cout << "Rendering complete. 100 pictures generated.";
	return 0;
}

void Update2(int index){
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;
	cout << "A 200pic movie will render in: " << 200 * (dt / 60000) << " minutes" << endl;
	if (index<100){
		lenses[0].focalLength = focals1[index];
		lenses[1].focalLength = focals2[index];
	}
	else if ((index >= 100) && (index<199)){
		float percent = ((float)index - 100.f) / 100.f;

		lenses[0].center = glm::vec3(glm::sin(2.f*PI*percent)*0.5, glm::cos(2.f*PI*percent)*0.5, -1.5f);
		lenses[1].center = lenses[0].center;
	}
	if (index == 100){
		lenses[0].focalLength = 2.0f;
		lenses[1].focalLength = 2.0f;
	}
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
	if (keystate[SDLK_i]) {
		lenseCenterStart -= down / float(4);
	}
	if (keystate[SDLK_j]) {
		lenseCenterStart -= right / float(4);
	}
	if (keystate[SDLK_k]) {
		lenseCenterStart += down / float(4);
	}
	if (keystate[SDLK_l]) {
		lenseCenterStart += right / float(4);
	}

	SetupLenses(lenseCenterStart);
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

/*
Initiates and places lenses at some point in the model
*/
void SetupLenses(vec3 center) {
	lenses = vector<Lense>(2);

	Lense lense;
	lense.center = center;
	lense.radius = 0.7f;
	lense.normal = glm::normalize(vec3(0, 0, -1));
	lense.focalLength = 4.f;
	lense.refractiveIndex = 1.5f;
	lense.index = 0;
	lenses[0] = lense;

	lense.center = center;
	lense.radius = 0.7f;
	lense.normal = glm::normalize(vec3(0, 0, 1));
	lense.focalLength = 3.f;
	lense.refractiveIndex = 1.5f;
	lense.index = 1;
	lenses[1] = lense;
}

float getSign(float number){
	if (number>0){
		return 1;
	}
	else{
		return -1;
	}
}

/*
Initiates noise of all lenses.
*/
void SetupLenseNoises(vector<Lense> lenses){
	lensenoises = vector<LenseNoiseMap>(2);
	float maxnoise = 9.8;
	float minimum = 0.2;
	float intmax = (float)INT_MAX;
	srand(time(NULL));
	for (int i = 0; i<lenses.size(); i++) {
		LenseNoiseMap noise;
		noise.maxangle = calculateMaxAngle(lenses[i]);
		noise.noisematrix = vector<vector<vec3>>(noisecells);

		for (int y = 0; y<noisecells; y++) {
			noise.noisematrix[y] = vector<vec3>(noisecells);
			for (int x = 0; x<noisecells; x++) {
				float rand1sign = getSign((float)(2 * (float)rand() / (float)intmax) - 1);
				float rand2sign = getSign((float)(2 * (float)rand() / (float)intmax) - 1);
				float rand3sign = getSign((float)(2 * (float)rand() / (float)intmax) - 1);
				float rand1 = (float)rand();
				float rand2 = (float)rand();
				float rand3 = (float)rand();
				noise.noisematrix[y][x] = vec3(((rand1 / intmax) + minimum)*maxnoise*rand1sign, ((rand2 / intmax) + minimum)*maxnoise*rand2sign, ((rand3 / intmax) + minimum)*maxnoise*rand3sign);
				//                cout << noise.noisematrix[y][x].x << "," << noise.noisematrix[y][x].y << "," << noise.noisematrix[y][x].z << "\n";
			}
			//            cout << "\n";
		}
		lensenoises[i] = noise;
	}
}

/*
Checks if a ray starting in point 'start' in direction 'dir' intersects a lense in its path.
Returns true of it does, otherwise false.
*/
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

			if (d1<0){
				d1 = d2;
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

/*
Calculates the refraction of the ray that is being traced. Calculates both the refraction that happens when the ray hits the lens
aswell as when it exits.
*/
void calculateRefraction(vec3 dirIn, vec3 lensePointIn, int lenseIndex, vec3& lensePointOut, vec3& dirOut) {
	vec3 insideRefractionDir;
	Lense lense = lenses[lenseIndex];
	calculateRefractionVector(lense, dirIn, lensePointIn, defaultRefractiveIndex, lense.refractiveIndex, insideRefractionDir);

	//The intersection when the ray exits the lense
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

/*
Calcuates a vector that describes the ray when it has hit the lense at the point 'pointIn'.
*/
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

	vec3 noise = calculateNoise(lense.index, glm::atan(normalVector.y / normalVector.z), glm::atan(normalVector.x / normalVector.z));

	vec3 rotatedVector = rotateVector(dirIn, pitch2 - pitch1, yaw1 - yaw2);

	//cout << "noise1: " << rotatedVector.x << "," << rotatedVector.y << "," << rotatedVector.z << "\n";

	dirOut = glm::normalize(rotatedVector) + noise;

	//    dirOut = dirIn;
	//cout << "noise2: " << dirOut.x << "," << dirOut.y << "," << dirOut.z << "\n";
}

/*
Rotates the 'vector' as specified by 'pith' and 'yaw'. Returnes the rotated vector.
*/
vec3 rotateVector(vec3 vector, float pitch, float yaw) {
	mat3 matPitch = mat3(vec3(1, 0, 0), vec3(0, cos(pitch), -sin(pitch)), vec3(0, sin(pitch), cos(pitch)));
	mat3 matYaw = mat3(vec3(cos(yaw), 0, sin(yaw)), vec3(0, 1, 0), vec3(-sin(yaw), 0, cos(yaw)));
	mat3 mat = matPitch * matYaw;
	vec3 rotated = mat * vector;
	return rotated;
}


float calculateMaxAngle(Lense lense){
	vec3 focalPoint = (-lense.focalLength)*lense.normal + lense.center;
	vec3 perpendicular;

	findPerpendicular(lense.normal, perpendicular);
	vec3 lensenormal = lense.normal;
	vec3 pointnormal = glm::normalize(perpendicular*lense.radius - focalPoint);

	return glm::acos(glm::dot(lensenormal, pointnormal) / (glm::length(lensenormal)*glm::length(pointnormal)));
}

vec3 calculateNoise(int lenseindex, float anglex, float angley){



	float distanceconstant = 1;
	float maxangle = lensenoises[lenseindex].maxangle;

	float interval = glm::abs(maxangle) * 2;
	float absanglex = maxangle + anglex;
	float absangley = maxangle + angley;

	float xquota = absanglex / interval;
	float yquota = absangley / interval;

	float xval = xquota*(noisecells - 1);
	float yval = yquota*(noisecells - 1);

	if (lenseindex == 0){
		//cout << xval << "," << yval << "\n";
		if ((xval>150) || (xval<0) || (yval>150) || (yval<0)){
			cout << "?????\n";
		}
	}
	vec3 noisesum = vec3(0, 0, 0);

	for (int x = 0; x<noisecells; x++){
		for (int y = 0; y<noisecells; y++) {



			float ydistance = glm::abs(y - yval);
			float xdistance = glm::abs(x - xval);
			float distance = glm::sqrt((xdistance*xdistance) + (ydistance*ydistance));

			float quota = distanceconstant / (distance*distance*distance*distance);
			noisesum = noisesum + (lensenoises[lenseindex].noisematrix[x][y] * quota);
		}
	}
	if (lenseindex == 1){
		//cout << yquota << "\n";
	}
	noisesum = glm::normalize(noisesum) / 100.f;
	//cout << "(" << noisesum.x << ", " << noisesum.y << ", " << noisesum.z << ")\n";

	return noisesum;
}

