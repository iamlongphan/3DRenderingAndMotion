#include "ofApp.h"
//Signed by Long Phan - 10/3/2023

// Intersect Ray with Plane  (wrapper on glm::intersect*
//
bool Plane::intersect(const Ray& ray, glm::vec3& point, glm::vec3&
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


// Convert (u, v) to (x, y, z) 
// We assume u,v is in [0, 1]
//
glm::vec3 ViewPlane::toWorld(float u, float v) {
	float w = width();
	float h = height();
	return (glm::vec3((u * w) + min.x, (v * h) + min.y, position.z));
}





// Get a ray from the current camera position to the (u, v) position on
// the ViewPlane
//
Ray RenderCam::getRay(float u, float v) {
	glm::vec3 pointOnPlane = view.toWorld(u, v);
	return(Ray(position, glm::normalize(pointOnPlane - position)));
}


//GetRaySamples for PointLights HERE.
int PointLight::getRaySamples(glm::vec3& p, vector<Ray>& raySampleHolder)
{
	//Light direction from the light that it is calling from
	glm::vec3 lightDirection = glm::normalize(this->position-p);

	//Here we're recalculating the origin, intersectionPoint + (FindingOriginalNormal and then moving that point slightly to avoid the shadow rounding errors), and then finally creating the ray.
	glm::vec3 shadowRayOrig = this->position;//+ (lightDirection * 0.1f);
	Ray shadowRay(shadowRayOrig, lightDirection);
	
	//We add it into the raySampleHolder so we can come back to this ray to later interpret the shadows and whether or not we add color in the main shading functions. (this will probably require a loop in the area light code)
	raySampleHolder.push_back(shadowRay);
	

	return raySampleHolder.size();
}
//END POINTLIGHT

//GetRaySamples for AREALIGHTS HERE
int AreaLight::getRaySamples(glm::vec3& p, vector<Ray>& raySampleHolder)
{
	//We gather the random points
	vector<glm::vec3> randomPoints = this->getRandomPoints();

	//From the random points, we then gather them to make rays from the point and then the plane, and we insert them to the vector
	for (glm::vec3 ranPt : randomPoints)
	{
		glm::vec3 lightDir = glm::normalize(ranPt-p);
		glm::vec3 shadowRayOrig = p + (lightDir - glm::normalize(p) * 0.1f);
		Ray shadowRay(shadowRayOrig, lightDir);
		raySampleHolder.push_back(shadowRay);
	}
	return raySampleHolder.size();
}
//END AREALIGHT

//Possible Errors That Were Solved (WEPTS): Dynamic casting objects caused a memory allocation error when trying to add spheres
ofColor ofApp::phong(glm::vec3& p, glm::vec3& norm, ofColor diff, ofColor spec, float pow)
{
	float lightintensity = intensitySlider;
	pow = powerslider;

	ofColor maxCol = ofColor::black;

	for (int z = 0; z < lights.size(); z++)
	{
		
		vector<Ray> rayVectorP;			
		int vectorSize = lights[z]->getRaySamples(p, rayVectorP);
		//Dont have to cast because its just gonna be one sample point if its just a regular light and then multiple points if it is an arealight.
		for (int i = 0; i < rayVectorP.size(); i++)
		{
			Ray shadowRay3 = rayVectorP.at(i);
			glm::vec3 lightDir = shadowRay3.d;
			bool shad3 = false;
			for (int j = 0; j < scene.size(); j++)
			{
				SceneObject* tempObj2 = scene[j];
				if (tempObj2->alight == false && tempObj2->intersect(shadowRay3, shadowRay3.p, glm::normalize(tempObj2->position)))
				{
					shad3 = true;
					break;
				}
			}
			if (!shad3)
			{
				float diffuseIntens = glm::max(glm::dot(norm, lightDir), 0.0f);
				//This is the color of our diffuseCoeff calculated from the light's color(white), diffuse of object color, diffuseIntensity * the lights intensity + slider intensity
				ofColor diffuseLight = (ofColor::white * diff * diffuseIntens * ((lights[z]->intensity * lightintensity) / glm::pow(glm::distance(p, lights[z]->position), 2))) / static_cast<float>(vectorSize);
				//This is the direction of our rendercam view
				glm::vec3 viewDir = glm::normalize(mainCam.getPosition());
				//Gives us the direction of the reflection.
				glm::vec3 reflect = glm::reflect(-lightDir, norm);
				//Intensity for specularLighting.
				float specularIntens = glm::pow(glm::max(glm::dot(viewDir, reflect), 0.0f), pow);
				//Specular color 
				ofColor specLight = (ofColor::white * spec * specularIntens * ((lights[z]->intensity * lightintensity) / glm::pow(glm::distance(p, lights[z]->position), 2))) / static_cast<float>(vectorSize); //removed  / static_cast<float>(vectorSize), testing on just diffuse color now.
				//Combine for Blinn-Phong
				maxCol += (diffuseLight + specLight) ;
			}
		}
	}
	return maxCol;
}

//WHY RAYTRACE COLORING 2? Its a band-aid fix for my motion blur portion of the if statements in the rayTraceColoring method, since we need to get the color multiple times to emulate jitter motion blur.
ofColor ofApp::rayTraceColoring2(Ray& r)
{
	ofColor colorForObj = ofColor::black;
	glm::vec3 normalAtInt;
	distance2 = std::numeric_limits<float>::infinity();
	SceneObject* closestObj = nullptr;
	float tempDist;
	glm::vec3 closestP;

	bool hit = false;
	for (int i = 0; i < scene.size(); i++)
	{
		glm::vec3 point;
		glm::vec3 normal;

		SceneObject* tempObj = scene[i];
		if (tempObj->intersect(r, point, normal))
		{
			hit = true;
			tempDist = glm::distance(r.p, point);
			if (tempDist < distance2)
			{

				distance2 = tempDist;
				closestObj = tempObj;
				normalAtInt = normal;
				closestP = point;

				if (hit)
				{
					if (closestObj->isfloor)
					{
						//This float is what changes it from either a 10x10 scale, or like a 50x50 etcetc
						float scale = 7.0f;

						//We need this in order to actually grab the height/width of the plane, because its currently a SceneObject*
						Plane* floorModel = dynamic_cast<Plane*>(closestObj);

						//Ucorrected and VCOrrected are found from respectively the X and Z of the floor plane and then scaled across the plane as necessary to tile. More explicitly, we take the x of the intersection point and the x from the vector of the floor plane
						float uCorrected = (closestP.x - floorModel->position.x) / floorModel->height;
						float vCorrected = (closestP.z - floorModel->position.z) / floorModel->width;
						uCorrected = uCorrected * scale;
						vCorrected = vCorrected * scale;

						//This is essentially our "wrapping" to ensure theres no aritfacts, we take just the float numbers so that floorX and likewise floorZ are getting correct values.
						float uTiledFloor = uCorrected - floor(uCorrected);

						//Clamp it to 0.0f and 1.0f to ensure it stays between these values, so my program doesn't crash.
						uTiledFloor = ofClamp(uTiledFloor, 0.0f, 1.0f);


						float floorX = fmod(uTiledFloor * (floorTextureImage.getHeight() - 0.5f), floorTextureImage.getHeight());

						float vTiledFloor = vCorrected - floor(vCorrected);
						vTiledFloor = ofClamp(vTiledFloor, 0.0f, 1.0f);


						float floorZ = fmod(vTiledFloor * (floorTextureImage.getWidth() - 0.5f), floorTextureImage.getWidth());

						ofColor floorDiff = floorTextureImage.getColor(floorX, floorZ);
						float floorXSpec = fmod(uTiledFloor * (floorTextureImageSPEC.getHeight() - 0.5f), floorTextureImageSPEC.getHeight());
						float floorZSpec = fmod(vTiledFloor * (floorTextureImageSPEC.getWidth() - 0.5f), floorTextureImageSPEC.getWidth());

						ofColor floorSpec = floorTextureImageSPEC.getColor(floorXSpec, floorZSpec);
						colorForObj = phong(closestP, normalAtInt, floorDiff, floorSpec, 100.0f);
					}
					else if (closestObj->iswall)
					{
						float scale = 6.0f;
						Plane* wallModel = dynamic_cast<Plane*>(closestObj);

						float uCorrected = (closestP.x - wallModel->position.x) / wallModel->height;
						float vCorrected = (closestP.y - wallModel->position.y) / wallModel->width;
						uCorrected = uCorrected * scale;
						vCorrected = vCorrected * scale;


						float uTiledWall = uCorrected - floor(uCorrected);

						uTiledWall = ofClamp(uTiledWall, 0.f, 1.f);

						float wallX = fmod(uTiledWall * (textureImage2.getHeight()) - 0.5f, textureImage2.getHeight());

						float vTiledWall = vCorrected - floor(vCorrected);
						vTiledWall = ofClamp(vTiledWall, 0.f, 1.f);

						float wallY = fmod(vTiledWall * (textureImage2.getWidth()) - 0.5f, textureImage2.getWidth());
						ofColor wallDiff = textureImage2.getColor(wallX, wallY);

						float wallXSpec = fmod(uTiledWall * (textureImageSPEC2.getHeight() - 0.5f), textureImageSPEC2.getHeight());
						float wallYSpec = fmod(vTiledWall * (textureImageSPEC2.getWidth() - 0.5f), textureImageSPEC2.getWidth());

						ofColor wallSpec = textureImageSPEC2.getColor(wallXSpec, wallYSpec);
						colorForObj = phong(closestP, normalAtInt, wallDiff, wallSpec, 100.0f);
					}
					//Sphere textureing for one sphere.
					else if (closestObj->isSphere1)
					{

						//Since we draw from its position and the radius expands anyway.. just use position 
						glm::vec3 sphPos = closestObj->position;
						glm::vec3 sphN = glm::normalize(closestP - sphPos);

						//angle from sphere normal's y in polar coordinates
						float upAng = glm::acos(-sphN.y);
						// Unit Circle
						float unitCircleAng = glm::atan(sphN.z, sphN.x);

						//Coordinates for sphere here.
						float u = unitCircleAng / (2.0f * PI) + 0.5f;
						float v = upAng / PI;

						ofColor sphDiff = sphereImage.getColor(u * sphereImage.getWidth(), v * sphereImage.getHeight());
						ofColor sphSpec = sphereImageSPEC.getColor(u * sphereImageSPEC.getWidth(), v * sphereImageSPEC.getHeight());

						colorForObj = phong(closestP, normalAtInt, sphDiff, sphSpec, 100.0f);
					}
					else if (closestObj->isSelected)
					{
						colorForObj = phong(closestP, normalAtInt, ofColor::white, closestObj->specularColor, 100.0f);
					}
					else if(!closestObj->alight)
					{
						colorForObj = phong(closestP, normalAtInt, closestObj->diffuseColor, closestObj->specularColor, 100.0f);
					}
				}

			}

		}
	}
	return colorForObj;
}



//RAY TRACING COLOR
ofColor ofApp::rayTraceColoring(Ray& r)
{
	ofColor colorForObj = ofColor::black;
	glm::vec3 normalAtInt;
	distance2 = std::numeric_limits<float>::infinity();
	SceneObject* closestObj = nullptr;
	float tempDist;
	glm::vec3 closestP;

	bool hit = false;
	for (int i = 0;  i < scene.size() ;i++)
	{
		glm::vec3 point;
		glm::vec3 normal;
		
		SceneObject* tempObj = scene[i];
		if (tempObj->intersect(r, point, normal))
		{
			hit = true;
			tempDist = glm::distance(r.p, point);
			if (tempDist < distance2)
			{
				distance2 = tempDist;
				closestObj = tempObj;
				normalAtInt = normal;
				closestP = point;
				
				if (hit)
				{
					//MOTION BLUR HERE 
					if (isAnimating && !closestObj->alight && closestObj == key1.obj2 && motionblurOnOff)
					{
						
						//Essentially what we're doing here is finding a lower and upper bound for the rx,ry,and rz, their purpose is to jitter the ray
						//In order to create the "blur" effect, now the reasoning for blur factor, is in the technique,
						//We're creating the min bound and max bound because as it jitters farther the warbling effect gets stronger,
						//so when we warble(jitter) closer, we get less blurr, and the effect is weaker.
						//This is important in especially the easeInterp which best uses the motion blur effect, because the linearInterp
						//Would be a constant velocity, therefore theres no real need to show a motion blur effect
						int middleOf2KeyFrames = (key2.frame + key1.frame)/2;
						//Blur fac here goes from a value of 0 to 1, but think of it as a bell curve with peak strength at the middle, thats the effect
						//That I wanted for blur factor, as it helps creates the bounds for the amount of blur I create.
						float blurFac = glm::clamp(1.0f - abs((float)currentFrame - middleOf2KeyFrames) / (middleOf2KeyFrames / 2), 0.0f, 1.0f);
						
						// makes the value loop back around  essentially, -0.0085 to 0.0085 are the bounds but we shrink and stretch as required
						//full equation of mix is  x * (1.0 - a) + y * a, where mix(x,y,a);
						float minB = glm::mix(-0.0001, -0.0085, blurFac);
						
						float maxB = glm::mix(0.0001, 0.0085, blurFac);
						
						//These are the values of which we jitter the direction of the ray to create the blurring effect, randomly in values within bounds.
						float rx = ofRandom(minB, maxB);
						float ry = ofRandom(minB, maxB);
						float rz = ofRandom(minB, maxB);
	
						//Create a jittery ray direction
						glm::vec3 jitteredRayDir = glm::normalize(glm::vec3(r.d.x + rx, r.d.y + ry, r.d.z + rz));
						Ray jitteredRay(r.p, jitteredRayDir);

						//We take the color of the jitter and reapply.
						colorForObj = rayTraceColoring2(jitteredRay);

					}
					//MOTION BLUR END
					//re-add else here for with motion blur
					 else if (closestObj->isfloor)
					{
						//This float is what changes it from either a 10x10 scale, or like a 50x50 etcetc
						float scale = 7.0f;

						//We need this in order to actually grab the height/width of the plane, because its currently a SceneObject*
						Plane* floorModel = dynamic_cast<Plane*>(closestObj);

						//Ucorrected and VCOrrected are found from respectively the X and Z of the floor plane and then scaled across the plane as necessary to tile. More explicitly, we take the x of the intersection point and the x from the vector of the floor plane
						float uCorrected = (closestP.x - floorModel->position.x) / floorModel->height;
						float vCorrected = (closestP.z - floorModel->position.z) / floorModel->width;
						uCorrected = uCorrected * scale;
						vCorrected = vCorrected * scale;

						//This is essentially our "wrapping" to ensure theres no aritfacts, we take just the float numbers so that floorX and likewise floorZ are getting correct values.
						float uTiledFloor = uCorrected - floor(uCorrected);

						//Clamp it to 0.0f and 1.0f to ensure it stays between these values, so my program doesn't crash.
						uTiledFloor = ofClamp(uTiledFloor, 0.0f, 1.0f);


						float floorX = fmod(uTiledFloor * (floorTextureImage.getHeight() - 0.5f), floorTextureImage.getHeight());

						float vTiledFloor = vCorrected - floor(vCorrected);
						vTiledFloor = ofClamp(vTiledFloor, 0.0f, 1.0f);


						float floorZ = fmod(vTiledFloor * (floorTextureImage.getWidth() - 0.5f), floorTextureImage.getWidth());

						ofColor floorDiff = floorTextureImage.getColor(floorX, floorZ);
						float floorXSpec = fmod(uTiledFloor * (floorTextureImageSPEC.getHeight() - 0.5f), floorTextureImageSPEC.getHeight());
						float floorZSpec = fmod(vTiledFloor * (floorTextureImageSPEC.getWidth() - 0.5f), floorTextureImageSPEC.getWidth());

						ofColor floorSpec = floorTextureImageSPEC.getColor(floorXSpec, floorZSpec);
						colorForObj = phong(closestP, normalAtInt, floorDiff, floorSpec, 100.0f);
					}
					else if (closestObj->iswall)
					{
						float scale = 6.0f;
						Plane* wallModel = dynamic_cast<Plane*>(closestObj);

						float uCorrected = (closestP.x - wallModel->position.x) / wallModel->height;
						float vCorrected = (closestP.y - wallModel->position.y) / wallModel->width;
						uCorrected = uCorrected * scale;
						vCorrected = vCorrected * scale;


						float uTiledWall = uCorrected - floor(uCorrected);

						uTiledWall = ofClamp(uTiledWall, 0.f, 1.f);

						float wallX = fmod(uTiledWall * (textureImage2.getHeight()) - 0.5f, textureImage2.getHeight());

						float vTiledWall = vCorrected - floor(vCorrected);
						vTiledWall = ofClamp(vTiledWall, 0.f, 1.f);

						float wallY = fmod(vTiledWall * (textureImage2.getWidth()) - 0.5f, textureImage2.getWidth());
						ofColor wallDiff = textureImage2.getColor(wallX, wallY);

						float wallXSpec = fmod(uTiledWall * (textureImageSPEC2.getHeight() - 0.5f), textureImageSPEC2.getHeight());
						float wallYSpec = fmod(vTiledWall * (textureImageSPEC2.getWidth() - 0.5f), textureImageSPEC2.getWidth());

						ofColor wallSpec = textureImageSPEC2.getColor(wallXSpec, wallYSpec);
						colorForObj = phong(closestP, normalAtInt, wallDiff, wallSpec, 100.0f);
					}
					//Selection code
					else if (closestObj->isSelected)
					{
						colorForObj = phong(closestP, normalAtInt, ofColor::white, closestObj->specularColor, 100.0f);
					}
					//Sphere Texturization //This is limited to any spheres that are tagged as such
					else if (closestObj->isSphere1)
					{
					
						//Since we draw from its position and the radius expands anyway.. just use position 
						glm::vec3 sphPos = closestObj->position;
						glm::vec3 sphN = glm::normalize(closestP - sphPos);

						//angle from sphere normal's y in polar coordinates
						float upAng = glm::acos(-sphN.y);
						// Unit Circle
						float theta = glm::atan(sphN.z, sphN.x);

						//Coordinates for sphere here.
						float u = theta / (2.0f * PI) + 0.5f;
						float v = upAng / PI;

						ofColor sphDiff = sphereImage.getColor(u * sphereImage.getWidth(), v * sphereImage.getHeight());
						ofColor sphSpec = sphereImageSPEC.getColor(u * sphereImageSPEC.getWidth(), v * sphereImageSPEC.getHeight());

						colorForObj = phong(closestP, normalAtInt, sphDiff, sphSpec, 100.0f);
					}
					else
					{
						colorForObj = phong(closestP, normalAtInt, closestObj->diffuseColor, closestObj->specularColor, 100.0f);
					}
				}
			}
		}
	}
	return colorForObj;
}

//END RAY TRACING COLOR

//Grab mainCam information to draw over scene and render
void ofApp::rayTrace() {
	//Camera grabs
	glm::vec3 cameraPos = mainCam.getPosition();
	glm::vec3 mainLook = mainCam.getLookAtDir();
	glm::vec3 mainUpVec = mainCam.getUpDir();
	float camFOV = mainCam.getFov();
	float camAR = mainCam.getAspectRatio();
	glm::mat4 viewMatrix = glm::lookAt(glm::vec3(cameraPos.x, cameraPos.y, cameraPos.z), glm::vec3(mainLook.x, mainLook.y, mainLook.z), glm::vec3(mainUpVec.x, mainUpVec.y, mainUpVec.z));
	glm::mat4 projMatrix = glm::perspective(glm::radians(camFOV), camAR, mainCam.getNearClip(), mainCam.getFarClip());
	ofColor color;
	for (int j = 0; j < imageHeight; j++) {
		for (int i = 0; i < imageWidth; i++) {

			float u = (i+0.5f) / imageWidth;
			float v = (j+0.5f) / imageHeight;	
			
			//NDC coord conversion
			float x = 2.0f * u - 1.0f;
			float y = 1.0f - 2.0f * v;

			//Nearfar Clip plane
			glm::vec4 clip(x, y, -1.0f, 1.0f);
			//Find where we are looking from from that * inverse of projection matrix
			glm::vec4 eye = glm::inverse(projMatrix) * clip;
			eye.z = -1.0f;
			eye.w = 0.0f;

			//Go back to world space to find the ray direction and then normalize for calculations.
			glm::vec3 rayD = glm::vec3(glm::inverse(viewMatrix) * eye);
			rayD = glm::normalize(rayD);
			Ray ray(cameraPos, rayD);
			color = rayTraceColoring(ray);
			image.setColor(i, j, color);
			
		}
	}

	
	image.update();
	//image.save("somerandomname.jpg");
}


//Adding spheres, self-explanatory.
void ofApp::addSphere()
{
	Sphere* newSphere = new Sphere(glm::vec3(-2, 1, -2), 1.f, ofColor::yellow);
	scene.push_back(newSphere);
}

//Deleting objects out of the scene or lights.
void ofApp::deleteObject(SceneObject* obj)
{
	
	if (obj->alight)
	{
		for (int i = 0; i < lights.size(); i++)
		{
			if (lights[i] == obj)
			{
				lights.erase(lights.begin() + i);
			}
		}
		for (int j = 0; j < scene.size(); j++)
		{
			if (scene[j] == obj)
			{
				scene.erase(scene.begin() + j);
			}
		}
	}
	else
	{
		for (int z = 0; z < scene.size(); z++)
		{
			if (scene[z] == obj)
			{
				scene.erase(scene.begin() + z);
			}
		}
	}

	
}
void ofApp::addLight()
{
	Light* newLight = new PointLight();
	newLight->position = glm::vec3(0, 6, 0);
	newLight->diffuseColor = ofColor::white;
	scene.push_back(newLight);
	lights.push_back(newLight);
}


//--------------------------------------------------------------
void ofApp::setup() {
	
	image.allocate(imageWidth, imageHeight, OF_IMAGE_COLOR);

	//We load all the textures in right here so we can just apply them asap.
	sphereImage.load("44_marble floor_DIFF.jpg");
	sphereImageSPEC.load("44_marble floor_SPEC.jpg");

	textureImage2.load("53_Floral wallpaper_DIFF.jpg");	
	textureImageSPEC2.load("53_Floral wallpaper_SPEC.jpg");

	floorTextureImage.load("32_Wood floor_DIFF.jpg");
	floorTextureImageSPEC.load("32_Wood floor_SPEC.jpg");
	
	ofSetBackgroundColor(ofColor::black);

	theCam = &mainCam;

	mainCam.setDistance(10);
	mainCam.setNearClip(.1);
	
	//Blue Sphere Now also the textured sphere.
	Sphere* a = new Sphere(glm::vec3(-3.f, 0.f, -6.5f), 1.5f, ofColor::blue);
	a->isSphere1 = true;
	//Green Sphere
	Sphere* b = new Sphere(glm::vec3(-1.f, -1.f, -3.f), 1.f, ofColor::darkGreen);
	b->isSphere1 = true;
	//Red Sphere
	Sphere* c = new Sphere(glm::vec3(1.f, 1.f, -7.f), 2.5f, ofColor::red);
	//Plane
	Plane* regPlane = new Plane(glm::vec3(0.f, -4.f, 0.f), glm::vec3(0, 1, 0), ofColor::cyan, 40, 40);
	regPlane->isfloor = true;
	Plane* regPlane2 = new Plane(glm::vec3(0, 0, -25.f), glm::vec3(0, 0, 1), ofColor::orange, 80, 80);
	regPlane2->iswall = true;
	
	//Light 1 PORPLE TOP
	PointLight* newLight = new PointLight();
	newLight->position = glm::vec3(2.f, 7.f, -2.f);
	newLight->intensity = 75;


	//Light 2 GREEN FRONTd
	PointLight* newLight2 = new PointLight();
	newLight2->position = glm::vec3(7.f, 4.f, 2.f);
	newLight2->intensity = 75;
	

	//Light 3 YELLOW LEFT
	PointLight* newLight3 = new PointLight();
	newLight3->position = glm::vec3(-6.f, 6.f, 1.f);
	newLight3->intensity = 50;

	//Push all objects into scene vector
	scene.push_back(c);
	scene.push_back(b);
	scene.push_back(a);
	scene.push_back(regPlane);
	scene.push_back(regPlane2);

	scene.push_back(newLight);
	scene.push_back(newLight2);
	scene.push_back(newLight3);
	lights.push_back(newLight);
	lights.push_back(newLight2);
	lights.push_back(newLight3);

	
	
	gui.setup();
	
	gui.add(motionblurOnOff.setup("Motion Blur", false));
	gui.add(intensitySlider.setup("Intensity Slider", 0.5, 0, 1));
	gui.add(powerslider.setup("Power Slider", 100, 10, 10000));
	
	gui.add(color_rSlider.setup("Red", 0, 0, 255));
	gui.add(color_gSlider.setup("Green", 0, 0, 255));
	gui.add(color_bSlider.setup("Blue", 0, 0, 255));
	gui.add(radiusController.setup("Radius", 2.f, 0, 10.f));
	

}

//--------------------------------------------------------------
void ofApp::update() {
	if (playBack) {
		nextFrame();
			
	}
	if (keyFramesSet() && (currentFrame >= key1.frame && currentFrame <= key2.frame))
	{
		//Animating flag important for motion blur
		isAnimating = true;
		//only using easeInterp because, motion blur doesn't really apply to linearInterp (constant velocity)
		glm::vec3 interpolatedPosition = easeInterp(currentFrame, key1.frame, key2.frame, key1.position, key2.position);
		key1.obj2->position = interpolatedPosition;
	}
	else
	{
		isAnimating = false;
	}

	if (radiusChange)
	{
		selected[0]->radius = radiusController;
	}

	rayTrace();
}

//--------------------------------------------------------------
void ofApp::draw() {

	theCam->begin();
	//For 3D depth.
	ofEnableDepthTest();

	ofNoFill();
	

	//Draw all scene objects.

	for (int i = 0; i < scene.size(); i++) 
	{
		if (objSelected() && scene[i] == selected[0])
		{
			if (!scene[i]->alight)
			{
				//Change color of object to RGB values in sliders.
				selected[0]->diffuseColor.b = color_bSlider;
				selected[0]->diffuseColor.r = color_rSlider;
				selected[0]->diffuseColor.g = color_gSlider;
			}
			//Selection flag
			scene[i]->isSelected = true;
		}
		else
		{
			scene[i]->isSelected = false;
		}
		scene[i]->draw();
	}
	

	ofDisableDepthTest();
	theCam->end();


	
	ofSetColor(ofColor::white);
	
	//Corner view vs Full view
	if (fullView)
		image.draw(0, 0, ofGetWidth(), ofGetHeight());
	else
		image.draw(ofGetWidth() - 300, 0, 300, 150);

	gui.draw();
	// output key frame info // Borrowed from Lab 4 
	string str1;
	str1 += "Frame: " + std::to_string(currentFrame) + " of " + std::to_string(frameMax - frameMin + 1);
	ofSetColor(ofColor::white);
	ofDrawBitmapString(str1, 5, 15);

}

//--------------------------------------------------------------
void ofApp::exit()
{
	//Destructor to save the ram
	for (SceneObject* object : scene)
	{
		delete object;
	}
	for (Light* light : lights)
	{
		delete light;
	}
}
void ofApp::keyPressed(int key) {
	switch (key) {
	case OF_KEY_F1:
		theCam = &mainCam;
		break;
	case OF_KEY_F2:
		break;
	case OF_KEY_F3:
		break;
	case OF_KEY_ALT:
		bAltKeyDown = true;
		if (!mainCam.getMouseInputEnabled())
			mainCam.enableMouseInput();
		break;
	case '.':
		nextFrame();
		break;
	case ',':
		prevFrame();
		break;
	case ' ':
		playBack = !playBack;
		break;
	case 'c':
		if (mainCam.getMouseInputEnabled())
			mainCam.disableMouseInput();
		else
			mainCam.enableMouseInput();
		break;
	case 'd':
		deleteObject(selected[0]);
		selected.clear();
		break;
	case 'e':
		resetKeyFrames();
		break;
	case 'f':
		fullView = !fullView;
		break;
	case 'g':
		//Screenshot button
		image.save("screenShot.jpg");
		break;
	case 'h':
		radiusChange = true;
		break;
	case 'i':
		lightIntensChange = true;
		break;
	case 'k':
		//Set keyframe 1/2
		if (objSelected()) setKeyFrame(selected[0]);
		break;
	
	case 'l':
		addLight();
		break;

	case'm':
		selected[0]->isSphere1 = true;
		break;
	
	case 'r':
		bScale = true;
		break;

	case 's':
		//The reason this WASN'T working was because it was a bunch of dynamic casts breaking stuff. never doing that again..
		addSphere();
		break;
	

	case 'x':
		rotateX = true;
		break;
	case 'y':
		rotateY = true;
		break;
	case 'z':
		rotateZ = true;
		break;
	default:
		break;
	}
	
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {
	switch (key) {
	case OF_KEY_ALT:
		bAltKeyDown = false;
		mainCam.disableMouseInput();
		break;
	case 'h':
		radiusChange = false;
		break;
	case 'i':
		lightIntensChange = false;
		break;
	case 'r':
		bScale = false;
		break;
	case 'x':
		rotateX = false;
		break;
	case 'y':
		rotateY = false;
		break;
	case 'z':
		rotateZ = false;
		break;
	default:
		break;
	}
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {
	//Borrowed from SkeletonBuilder
	if (objSelected() && bDrag)
	{
		glm::vec3 point;
		mouseToDragPlane(x, y, point);
		if (rotateX)
		{
			selected[0]->rotation += glm::vec3((point.x - lastPoint.x) * 20.0, 0, 0);

		}
		else if (rotateY)
		{
			selected[0]->rotation += glm::vec3(0, (point.x - lastPoint.x) * 20.0, 0);
		}
		else if (rotateZ)
		{
			selected[0]->rotation += glm::vec3(0, 0, (point.x - lastPoint.x) * 20.0);
		}
		else if (bScale) {
			float scalef = (point.x - lastPoint.x) * .1;
			selected[0]->scale.y += scalef;
			selected[0]->scale.x += scalef / 2;
			selected[0]->scale.z += scalef / 2;
		}
		else
		{
			selected[0]->position += (point - lastPoint);
		}
		lastPoint = point;
	}
}

//--------------------------------------------------------------
//Code Borrowed From SkeletonBuilder for selection.
void ofApp::mousePressed(int x, int y, int button) {
	if(mainCam.getMouseInputEnabled()) 
		return;
	
	selected.clear();

	//
	// test if something selected
	//
	vector<SceneObject*> hits;

	glm::vec3 p = theCam->screenToWorld(glm::vec3(x, y, 0));
	glm::vec3 d = p - theCam->getPosition();
	glm::vec3 dn = glm::normalize(d);

	// check for selection of scene objects
	//
	for (int i = 0; i < scene.size(); i++) {

		glm::vec3 point, norm;

		//  We hit an object
		//
		if (scene[i]->isSelectable && scene[i]->intersect(Ray(p, dn), point, norm)) {
			hits.push_back(scene[i]);
		}
	}


	// if we selected more than one, pick nearest
	//
	SceneObject* selectedObj = NULL;
	if (hits.size() > 0) {
		selectedObj = hits[0];
		float nearestDist = std::numeric_limits<float>::infinity();
		for (int n = 0; n < hits.size(); n++) {
			float dist = glm::length(hits[n]->position - theCam->getPosition());
			if (dist < nearestDist) {
				nearestDist = dist;
				selectedObj = hits[n];
			}
		}
	}
	if (selectedObj) {
		selected.push_back(selectedObj);
		bDrag = true;
		mouseToDragPlane(x, y, lastPoint);

	}
	else {
		selected.clear();
	}
}

//  This projects the mouse point in screen space (x, y) to a 3D point on a plane
//  normal to the view axis of the camera passing through the point of the selected object.
//  If no object selected, the plane passing through the world origin is used.
//	BORROWED FROM SKELETON BUILDER 
bool ofApp::mouseToDragPlane(int x, int y, glm::vec3& point) {
	glm::vec3 p = theCam->screenToWorld(glm::vec3(x, y, 0));
	glm::vec3 d = p - theCam->getPosition();
	glm::vec3 dn = glm::normalize(d);

	float dist;
	glm::vec3 pos;
	if (objSelected()) {
		pos = selected[0]->position;
	}
	else pos = glm::vec3(0, 0, 0);
	if (glm::intersectRayPlane(p, dn, pos, glm::normalize(theCam->getZAxis()), dist)) {
		point = p + dn * dist;
		return true;
	}
	return false;
}
//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {
	bDrag = false;
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}
