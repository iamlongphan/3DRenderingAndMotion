//
//  RayCaster - Set of simple classes to create a camera/view setup for our Ray Tracer HW Project
//
//  I've included these classes as a mini-framework for our introductory ray tracer.
//  You are free to modify/change.   
//
//  These classes provide a simple render camera which can can return a ray starting from
//  it's position to a (u, v) coordinate on the view plane.
//
//  The view plane is where we can locate our photorealistic image we are rendering.
//  The field-of-view of the camera by moving it closer/further 
//  from the view plane.  The viewplane can be also resized.  When ray tracing an image, the aspect
//  ratio of the view plane should the be same as your image. So for example, the current view plane
//  default size is ( 6.0 width by 4.0 height ).   A 1200x800 pixel image would have the same
//  aspect ratio.
//
//  This is not a complete ray tracer - just a set of skelton classes to start.  The current
//  base scene object only stores a value for the diffuse/specular color of the object (defaut is gray).
//  at some point, we will want to replace this with a Material class that contains these (and other 
//  parameters)
//  
//  (c) Kevin M. Smith  - 24 September 2018
//


//Signed by Long Phan / 10/3/2023
#pragma once

#include "ofMain.h"
#include <glm/gtx/intersect.hpp>
#include <ofxGui.h>
//#include <ofxOpenCv.h>
#include "glm/gtx/euler_angles.hpp"





//  General Purpose Ray class 
//
class Ray {
public:
	Ray(glm::vec3 p, glm::vec3 d) { this->p = p; this->d = d; }
	void draw(float t) { ofDrawLine(p, p + t * d); }
	
	glm::vec3 evalPoint(float t) {
		return (p + t * d);
	}

	

	glm::vec3 p, d;
};


//  Base class for any renderable object in the scene
//
class SceneObject {
public: 
	virtual void draw() = 0;    // pure virtual funcs - must be overloaded
	virtual bool intersect(const Ray &ray, glm::vec3 &point, glm::vec3 &normal) { cout << "SceneObject::intersect" << endl; return false; }
	
	// commonly used transformations
//
	glm::mat4 getRotateMatrix() {
		return (glm::eulerAngleYXZ(glm::radians(rotation.y), glm::radians(rotation.x), glm::radians(rotation.z)));   // yaw, pitch, roll 
	}
	glm::mat4 getTranslateMatrix() {
		return (glm::translate(glm::mat4(1.0), glm::vec3(position.x, position.y, position.z)));
	}
	glm::mat4 getScaleMatrix() {
		return (glm::scale(glm::mat4(1.0), glm::vec3(scale.x, scale.y, scale.z)));
	}


	glm::mat4 getLocalMatrix() {

		// get the local transformations + pivot
		//
		glm::mat4 scale = getScaleMatrix();
		glm::mat4 rotate = getRotateMatrix();
		glm::mat4 trans = getTranslateMatrix();

		// handle pivot point  (rotate around a point that is not the object's center)
		//
		glm::mat4 pre = glm::translate(glm::mat4(1.0), glm::vec3(-pivot.x, -pivot.y, -pivot.z));
		glm::mat4 post = glm::translate(glm::mat4(1.0), glm::vec3(pivot.x, pivot.y, pivot.z));



		return (trans * post * rotate * pre * scale);

	}

	glm::mat4 getMatrix() {

		// if we have a parent (we are not the root),
		// concatenate parent's transform (this is recursive)
		// 
		if (parent) {
			glm::mat4 M = parent->getMatrix();
			return (M * getLocalMatrix());
		}
		else return getLocalMatrix();  // priority order is SRT
	}

	// get current Position in World Space
	//
	glm::vec3 getPosition() {
		return (getMatrix() * glm::vec4(0.0, 0.0, 0.0, 1.0));
	}

	// set position (pos is in world space)
	//
	void setPosition(glm::vec3 pos) {
		position = glm::inverse(getMatrix()) * glm::vec4(pos, 1.0);
	}

	// return a rotation  matrix that rotates one vector to another
	//
	glm::mat4 rotateToVector(glm::vec3 v1, glm::vec3 v2);

	SceneObject* parent = NULL;        // if parent = NULL, then this obj is the ROOT
	vector<SceneObject*> childList;

	// rotate pivot
	//
	glm::vec3 pivot = glm::vec3(0, 0, 0);

	// any data common to all scene objects goes here
	glm::vec3 position = glm::vec3(0, 0, 0);
	glm::vec3 rotation = glm::vec3(0, 0, 0);
	glm::vec3 scale = glm::vec3(1, 1, 1);
	float radius = 1.0f;
	//glm::vec3 posAt(float time);
	// material properties (we will ultimately replace this with a Material class - TBD)
	//
	ofColor diffuseColor = ofColor::white;    // default colors - can be changed.
	ofColor specularColor = ofColor::lightGray;
	//Figure if its a light oo not
	bool alight = false;
	bool aplane = false;
	bool iswall = false;
	bool isfloor = false;
	bool isSelectable = true;
	bool isSelected = false;
	bool isSphere1 = false;
};

class KeyFrame2 {
public:
	int frame = -1;     //  -1 => no key is set;
	glm::vec3 position = glm::vec3(0, 0, 0);   // translate channel
	glm::vec3 rotation = glm::vec3(0, 0, 0);   // rotate channel
	glm::vec3 scale = glm::vec3(1, 1, 1);   // rotate channel
	SceneObject* obj2 = NULL;               // object that is keyframed
};


class Light : public SceneObject {
public:
	float intensity = 0.0;
	bool isSelectable = true;
	void draw()
	{
		ofSetColor(this->diffuseColor);
		ofDrawSphere(position, radius);
	}
	
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, radius, point, normal));
	}
	virtual int getRaySamples(glm::vec3& p, vector<Ray>& sampleHolder) = 0;
};



//PointLight class here
class PointLight : public Light {
public: 

	PointLight()
	{
		intensity = 25.f;
		alight = true;
		isSelectable = true;
		radius = 0.5f;
	}
	//Overloaded this in .cpp
	int getRaySamples(glm::vec3& p, vector<Ray>& sampleHolder);

	void draw()
	{
		ofDrawSphere(position, .5);
	}
};

class AreaLight : public Light {

public:
	float width = 5, height = 5; //required
	//glm::vec3 position = glm::vec3(0,0,0); //required 
	int nDivsWidth = 10, nDivsHeight = 10; // required
	int nSamples = 1; //required

	glm::vec3 up = glm::vec3(0, 0, 1); //Axis for Z
	glm::vec3 positiveRight = glm::vec3(1, 0, 0); //Axis for X
	glm::vec3 normal = glm::vec3(0, -1, 0);
	
	//Default Constructor
	AreaLight()
	{
		float width = 5;
		float height = 5;
		position = glm::vec3(0, 6, 0);
		int nDivsWidth = 10;
		int nDivsHeight = 10;
		int nSamples = 1;
		alight = true;
		isSelectable = true;
		
	}

	//Constructor for if we want to be speciic
	AreaLight(float w, float h, glm::vec3 pos, int nDivsW, int nDivsH, int samples)
	{
		this->width = w;
		this->height = h;
		this->position = pos;
		this->nDivsWidth = nDivsW;
		this->nDivsHeight = nDivsH;
		this->nSamples = samples;
		alight = true;
		isSelectable = true;
	}


	int getRaySamples(glm::vec3& p, vector<Ray>& rayVector);


	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3&
		normalAtIntersect) {
		float dist;
		bool insidePlane = false;
		bool hit = glm::intersectRayPlane(ray.p, ray.d, position, this->normal,
			dist);
		if (hit) {
			Ray r = ray;
			point = r.evalPoint(dist);
			normalAtIntersect = this->normal;
			glm::vec2 xrange = glm::vec2(position.x - width / 2, position.x + width
				/ 2);
			glm::vec2 yrange = glm::vec2(position.y - width / 2, position.y + width
				/ 2);
			glm::vec2 zrange = glm::vec2(position.z - height / 2, position.z +
				height / 2);
			// horizontal
			//
			if (normal == glm::vec3(0, 1, 0) || normal == glm::vec3(0, -1, 0)) {
				if (point.x < xrange[1] && point.x > xrange[0] && point.z <
					zrange[1] && point.z > zrange[0]) {
					insidePlane = true;
				}
			}
			// front or back
			//
			else if (normal == glm::vec3(0, 0, 1) || normal == glm::vec3(0, 0, -1))
			{
				if (point.x < xrange[1] && point.x > xrange[0] && point.y <
					yrange[1] && point.y > yrange[0]) {
					insidePlane = true;
				}
			}
			// left or right
			//
			else if (normal == glm::vec3(1, 0, 0) || normal == glm::vec3(-1, 0, 0))
			{
				if (point.y < yrange[1] && point.y > yrange[0] && point.z <
					zrange[1] && point.z > zrange[0]) {
					insidePlane = true;
				}
			}
		}
		return insidePlane;
	}

	//Testing function for getting random points
	glm::vec3 getRandomPoint()
	{
		float wid = ofRandomuf();
		float hei = ofRandomuf();
		return this->position + this->positiveRight * (-0.5f + wid) * this->width + this->up * (-0.5f + hei) * this->height;
	}

	vector<glm::vec3> getRandomPoints()
	{
		vector<glm::vec3> randomPoints;
		//Because we're looking at a grid, we grab the reciprocal in order to also get the negative ends of a grid (think of a coordinate plane... a plane 5x5 would be 2.5 in the x and -2.5, and 2.5 in the z and -2.5 in the z as well so we must consider that.
		float invXnDivs = 1.0 / this->nDivsWidth;
		float invZnDivs = 1.0 / this->nDivsHeight;
		//Essentially we cell by cell from the x axis first.
		for (int row = 0; row < this->nDivsHeight; row++)
		{
			for (int col = 0; col < this->nDivsWidth; col++)
			{
				for (int samples = 0; samples < nSamples; samples++)
				{
					//These return a random float between 0-1, these floats are made just to randomize the points that we're gettinng
					float wid = ofRandomuf();
					float hei = ofRandomuf();
					// The calculations are made so that instead of getting all the points from maybe the edge of the cell or the center of the cell we distribute them across the cells that we access and as such thats what we have.
					glm::vec3 pointOnRegion = this->position + this->positiveRight * ((col + 0.5f) * invXnDivs + (wid - 0.5f) * invXnDivs) * this->width + this->up * ((row + 0.5f) * invZnDivs + (hei - 0.5f) * invZnDivs - 0.5f) * this->height; 
					//cout << pointOnRegion << endl;
					randomPoints.push_back(pointOnRegion);
				}
			}
		}
		return randomPoints;
	}
	
	void draw()
	{
		ofPlanePrimitive plane2;
		//ofSetColor(this->diffuseColor);
		plane2.setPosition(this->position);
		plane2.setWidth(this->width);
		plane2.setHeight(this->height);
		plane2.setResolution(this->nDivsWidth+1, this->nDivsHeight+1);
		plane2.rotateDeg(90, 1, 0, 0);
		plane2.drawWireframe();
		glm::vec3 normal(0.f, -1.f, 0.f);
	}


};

//  General purpose sphere  (assume parametric)
//
class Sphere: public SceneObject {
public:
	Sphere(glm::vec3 p, float r, ofColor diffuse = ofColor::lightGray) { position = p; radius = r; diffuseColor = diffuse; }
	Sphere() {}
	bool intersect(const Ray &ray, glm::vec3 &point, glm::vec3 &normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, radius, point, normal));
	}
	
	void draw()  { 
		//ofSetColor(diffuseColor);
		
		ofDrawSphere(position, radius); 
	}
	

};

//  Mesh class (will complete later- this will be a refinement of Mesh from Project 1)
//
class Mesh : public SceneObject {
	bool intersect(const Ray &ray, glm::vec3 &point, glm::vec3 &normal) { return false;  }
	void draw() { }
	
};


//  General purpose plane 
//
class Plane : public SceneObject {
public:
	Plane(glm::vec3 p, glm::vec3 n, ofColor diffuse = ofColor::green, float w =
		40, float h = 40) {
		position = p; normal = n;
		width = w;
		height = h;
		diffuseColor = diffuse;
		aplane = true;
		isSelectable = false;
		if (normal == glm::vec3(0, 1, 0))
			plane.rotateDeg(-90, 1, 0, 0);
		else if (normal == glm::vec3(0, -1, 0))
			plane.rotateDeg(90, 1, 0, 0);
		else if (normal == glm::vec3(1, 0, 0))
			plane.rotateDeg(90, 0, 1, 0);
		else if (normal == glm::vec3(-1, 0, 0))
			plane.rotateDeg(-90, 0, 1, 0);
	}
	Plane() {
		normal = glm::vec3(0, 1, 0);
		plane.rotateDeg(90, 1, 0, 0);
		isSelectable = false;
		aplane = true;
	}
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal);
	float sdf(const glm::vec3& p);
	glm::vec3 getNormal(const glm::vec3& p) { return this->normal; }
	void draw() {
		plane.setPosition(position);
		plane.setWidth(width);
		plane.setHeight(height);
		plane.setResolution(4, 4);
		//plane.drawWireframe(); 
		ofSetColor(this->diffuseColor);
		plane.draw();
	}
	ofPlanePrimitive plane;
	glm::vec3 normal;
	float width = 20;
	float height = 20;
};


// view plane for render camera
// 
class  ViewPlane: public Plane {
public:
	ViewPlane(glm::vec2 p0, glm::vec2 p1) { min = p0; max = p1; }

	ViewPlane() {                         // create reasonable defaults (6x4 aspect)
		min = glm::vec2(-2, -1);
		max = glm::vec2(2, 1);			
		position = glm::vec3(0, 0, 5);
		normal = glm::vec3(0, 0, 1);      // viewplane currently limited to Z axis orientation
	}

	void setSize(glm::vec2 min, glm::vec2 max) { this->min = min; this->max = max; }
	float getAspect() { return width() / height(); }

	glm::vec3 toWorld(float u, float v);   //   (u, v) --> (x, y, z) [ world space ]

	void draw() {
		ofNoFill();
		ofDrawRectangle(glm::vec3(min.x, min.y, position.z), width(), height());
	}
	

	
	float width() {
		return (max.x - min.x);
	}
	float height() {
		return (max.y - min.y); 
	}

	// some convenience methods for returning the corners
	//
	glm::vec2 topLeft() { return glm::vec2(min.x, max.y); }
	glm::vec2 topRight() { return max; }
	glm::vec2 bottomLeft() { return min;  }
	glm::vec2 bottomRight() { return glm::vec2(max.x, min.y); }

	//  To define an infinite plane, we just need a point and normal.
	//  The ViewPlane is a finite plane so we need to define the boundaries.
	//  We will define this in terms of min, max  in 2D.  
	//  (in local 2D space of the plane)
	//  ultimately, will want to locate the ViewPlane with RenderCam anywhere
	//  in the scene, so it is easier to define the View rectangle in a local'
	//  coordinate system.
	//
	glm::vec2 min, max;
};


//  render camera  - currently must be z axis aligned (we will improve this in project 4)
//
class RenderCam: public SceneObject {
public:
	RenderCam() {
		position = glm::vec3(0, 0, 5);
		aim = glm::vec3(0, 0, -1);
	}
	Ray getRay(float u, float v);
	void draw() { ofDrawBox(position, 1.0); };
	
	void drawFrustum();

	glm::vec3 aim;
	ViewPlane view;          // The camera viewplane, this is the view that we will render 
};

  

class ofApp : public ofBaseApp {

public:
	void setup();
	void update();
	void draw();
	void exit();
	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void mouseEntered(int x, int y);
	void mouseExited(int x, int y);
	void windowResized(int w, int h);
	void dragEvent(ofDragInfo dragInfo);
	void gotMessage(ofMessage msg);
	//void drawGrid();
	//void drawAxis(glm::vec3 position);


	//RayTrace/Render Methods
	void rayTrace();
	ofColor rayTraceColoring(Ray& r);
	ofColor rayTraceColoring2(Ray& r);
	ofColor phong(glm::vec3& p, glm::vec3& norm, ofColor diffuse, ofColor specular, float power);
	//ofColor lambert(glm::vec3& p, glm::vec3& norm, ofColor diffuse);
	float distance2;
	vector<SceneObject*> scene;
	vector<Light*> lights;
	

	// Show or Hide Images
	bool bHide = true;
	bool bShowImage = false;


	//Camera Related Stuff
	ofEasyCam  mainCam;
	ofCamera* theCam;    // set to current camera either mainCam or sideCam


	// Image Related Stuff - Mainly Render scene + Textures
	ofImage image;
	ofImage sphereImage;
	ofImage sphereImageSPEC;
	ofImage textureImage2;
	ofImage textureImageSPEC2;
	ofImage floorTextureImage;
	ofImage floorTextureImageSPEC;
	bool fullView = true;



	

	// ofxGui Related Stuff //Probably will want to change the power slider and intensity slider to individual things for the lights and then we'll deal with the colors
	ofxPanel gui;
	ofxFloatSlider powerslider;
	ofxFloatSlider intensitySlider;
	ofxToggle motionblurOnOff;
	ofxIntSlider color_rSlider;
	ofxIntSlider color_gSlider;
	ofxIntSlider color_bSlider;

	ofxFloatSlider radiusController;
	bool radiusChange;

	
	// Mouse Hit Ray Related Stuff
	glm::vec3 mouseDir;
	glm::vec3 mouseOrigin;
	glm::vec3 mouseLast;
	glm::vec3 intersect;
	glm::vec3 hitPoint;
	float distance;
	bool hitPlane = false;
	int imageWidth = 300;
	int imageHeight = 200;
	bool PlaneIntersect = false;


	// KeyFrame Animation Related Stuff
	bool playBack = false;
	bool bKey2Next = false;
	bool bScale = false;
	int currentFrame = 1;
	int frameMax = 100;
	int frameMin = 1;
	KeyFrame2 key1, key2;
	bool isAnimating = false;

	
	void nextFrame() {
		currentFrame = (currentFrame == frameMax ? frameMin : currentFrame + 1);
	}
	void prevFrame() {
		currentFrame = (currentFrame == frameMin ? currentFrame : currentFrame - 1);
	}
	void startPlayback() {
		playBack = true;
	}

	void stopPlayback() {
		playBack = false;
	}

	//// check if both keys set
	bool keyFramesSet() { return (key1.frame != -1 && key2.frame != -1); }

	//// linear interpolation between two keyframes
	glm::vec3 linearInterp(int frame, int frameStart, int frameEnd, const glm::vec3& start, const glm::vec3& end) {
		return mapVec(frame, frameStart, frameEnd, start, end);
	}

	//// example non-linear interpolation between two keyframes (ease-in ease-out)
	////
	glm::vec3 easeInterp(int frame, int frameStart, int frameEnd, const glm::vec3& start, const glm::vec3& end) {

		// normalize range (0 to 1) and input to ease formula
	//	//
		float s = ease(ofMap(frame, frameStart, frameEnd, 0.0, 1.0));

		return mapVec(s, 0.0, 1.0, start, end);
	}


	////  ease-in and ease-out interpolation between two key frames
	////  this function produces a sigmoid curve normalized in x, y in (0 to 1);
	////
	float ease(float x) {
		return (x * x / (x * x + (1 - x) * (1 - x)));
	}

	//// helper functions to use ofMap on a vector
	////
	//// input a float value in a float range, output a vector 
	////
	glm::vec3 mapVec(float val, float start, float end, const glm::vec3& outStart, const glm::vec3& outEnd) {
		return glm::vec3(
			ofMap(val, start, end, outStart.x, outEnd.x),
			ofMap(val, start, end, outStart.y, outEnd.y),
			ofMap(val, start, end, outStart.z, outEnd.z));
	}

	// input a vec3 value in a vec3 range, output a vector 
	//
	glm::vec3 mapVec(const glm::vec3& val, const glm::vec3& start, const glm::vec3& end, const glm::vec3& outStart, const glm::vec3& outEnd) {
		return glm::vec3(
			ofMap(val.x, start.x, end.x, outStart.x, outEnd.x),
			ofMap(val.y, start.y, end.y, outStart.y, outEnd.y),
			ofMap(val.z, start.z, end.z, outStart.z, outEnd.z));
	}

	// set keyframe for SceneObject at current frame
	// call this function the first time and key1 is set
	// call this function again and key2 is set.
	// this "cycles" until you call resetKeyFrames();
	//
	void setKeyFrame(SceneObject* obj) {
		if (bKey2Next && key1.obj2 == obj) {
			key2.frame = currentFrame;
			key2.obj2 = obj;
			key2.position = obj->position;
			key2.rotation = obj->rotation;
			key2.scale = obj->scale;
			bKey2Next = false;
		}
		else {
			key1.frame = currentFrame;
			key1.obj2 = obj;
			key1.position = obj->position;
			key1.rotation = obj->rotation;
			key1.scale = obj->scale;
			bKey2Next = true;
		}
	}

	// reset key frames
	//
	void resetKeyFrames() {
		key1.frame = key2.frame = -1;
		key1.obj2 = key2.obj2 = NULL;
		bKey2Next = false;
	}



	// Object Related Stuff
	vector<SceneObject*> selected;
	bool mouseToDragPlane(int x, int y, glm::vec3& point);
	bool objSelected() { return (selected.size() ? true : false); };
	void addSphere();
	void deleteObject(SceneObject* obj);
	void addLight();
	bool bDrag = false;
	bool bAltKeyDown = false;
	bool rotateX = false;
	bool rotateY = false;
	bool rotateZ = false;
	glm::vec3 lastPoint;
	bool lightIntensChange = false;
	
};
 