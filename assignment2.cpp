#include "stdafx.h"

#include <iostream>
#include <windows.h>
#include <gl.h>
#include <glu.h>
#include <glut.h>
#include <math.h>
#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <stdio.h>										// Standard I/O header 
#include <stdlib.h>

using namespace std;

#define PI			3.1415926
#define	COLORNUM		15
float	ColorArr[COLORNUM][3] = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 },
{ 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 },
{ 0.3, 0.3, 0.3 }, { 0.5, 0.5, 0.5 }, { 0.9, 0.9, 0.9 },
{ 1.0, 0.5, 0.5 }, { 0.5, 1.0, 0.5 }, { 0.5, 0.5, 1.0 },
{ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 128.0 / 255, 192.0 / 255 } };
#define PI		3.1415926
int	angle = 0;
bool check = false;
bool check_arm1 = false;
bool check_arm2 = false;
bool check_scroll = false;
bool check_r = false;
bool check_g = false;
bool check_b = false;
float arm1_r = 0.0, arm1_g = 1.0, arm1_b = 1.0;
float arm2_r = 1.0, arm2_g = 0.0, arm2_b = 1.0;
float scroll_r = 1.0, scroll_g = 1.0, scroll_b = 0.0;
bool check_shadow1 = false;
bool check_shadow2 = false;
float r, g, b;
#define __TGA_H__




typedef	struct
{
	GLubyte	* imageData;									// Image Data (Up To 32 Bits)
	GLuint	bpp;											// Image Color Depth In Bits Per Pixel
	GLuint	width;											// Image Width
	GLuint	height;											// Image Height
	GLuint	texID;											// Texture ID Used To Select A Texture
	GLuint	type;											// Image Type (GL_RGB, GL_RGBA)
} Texture;

typedef struct
{
	GLubyte Header[12];									// TGA File Header
} TGAHeader;


typedef struct
{
	GLubyte		header[6];								// First 6 Useful Bytes From The Header
	GLuint		bytesPerPixel;							// Holds Number Of Bytes Per Pixel Used In The TGA File
	GLuint		imageSize;								// Used To Store The Image Size When Setting Aside Ram
	GLuint		temp;									// Temporary Variable
	GLuint		type;
	GLuint		Height;									//Height of Image
	GLuint		Width;									//Width ofImage
	GLuint		Bpp;									// Bits Per Pixel
} TGA;

Texture floorTex;
TGAHeader tgaheader;									// TGA header
TGA tga;										// TGA image data



GLubyte uTGAcompare[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };	// Uncompressed TGA Header


bool LoadTGA(Texture * texture, char * filename)				// Load a TGA file
{
	FILE * fTGA;												// File pointer to texture file
	fopen_s(&fTGA, filename, "rb");								// Open file for reading

	if (fTGA == NULL)											// If it didn't open....
	{
		return false;														// Exit function
	}

	if (fread(&tgaheader, sizeof(TGAHeader), 1, fTGA) == 0)					// Attempt to read 12 byte header from file
	{
		if (fTGA != NULL)													// Check to seeiffile is still open
		{
			fclose(fTGA);													// If it is, close it
		}
		return false;														// Exit function
	}

	// an Uncompressed TGA image
	if (fread(tga.header, sizeof(tga.header), 1, fTGA) == 0)					// Read TGA header
	{
		if (fTGA != NULL)													// if file is still open
		{
			fclose(fTGA);													// Close it
		}
		return false;														// Return failular
	}

	texture->width = tga.header[1] * 256 + tga.header[0];					// Determine The TGA Width	(highbyte*256+lowbyte)
	texture->height = tga.header[3] * 256 + tga.header[2];					// Determine The TGA Height	(highbyte*256+lowbyte)
	texture->bpp = tga.header[4];										// Determine the bits per pixel
	tga.Width = texture->width;										// Copy width into local structure						
	tga.Height = texture->height;										// Copy height into local structure
	tga.Bpp = texture->bpp;											// Copy BPP into local structure

	if ((texture->width <= 0) || (texture->height <= 0) || ((texture->bpp != 24) && (texture->bpp != 32)))	// Make sure all information is valid
	{
		if (fTGA != NULL)													// Check if file is still open
		{
			fclose(fTGA);													// If so, close it
		}
		return false;														// Return failed
	}

	if (texture->bpp == 24)													// If the BPP of the image is 24...
		texture->type = GL_RGB;											// Set Image type to GL_RGB
	else																	// Else if its 32 BPP
		texture->type = GL_RGBA;											// Set image type to GL_RGBA

	tga.bytesPerPixel = (tga.Bpp / 8);									// Compute the number of BYTES per pixel
	tga.imageSize = (tga.bytesPerPixel * tga.Width * tga.Height);		// Compute the total amout ofmemory needed to store data
	texture->imageData = (GLubyte *)malloc(tga.imageSize);					// Allocate that much memory

	if (texture->imageData == NULL)											// If no space was allocated
	{
		fclose(fTGA);														// Close the file
		return false;														// Return failed
	}

	if (fread(texture->imageData, 1, tga.imageSize, fTGA) != tga.imageSize)	// Attempt to read image data
	{
		if (texture->imageData != NULL)										// If imagedata has data in it
		{
			free(texture->imageData);										// Delete data from memory
		}
		fclose(fTGA);														// Close file
		return false;														// Return failed
	}

	// switch R and B
	for (int i = 0; i < tga.imageSize; i += tga.bytesPerPixel)
	{
		GLubyte temp = texture->imageData[i];
		texture->imageData[i] = texture->imageData[i + 2];
		texture->imageData[i + 2] = temp;
	}


	fclose(fTGA);															// Close file
	return true;															// All went well, continue on
}

void loadTextures(void)	{
	bool status = LoadTGA(&floorTex, "MARBLE.tga");
	if (status) {
		glGenTextures(1, &floorTex.texID);
		glBindTexture(GL_TEXTURE_2D, floorTex.texID);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, floorTex.width,
			floorTex.height, 0,
			GL_RGB, GL_UNSIGNED_BYTE, floorTex.imageData);

		if (floorTex.imageData)
			free(floorTex.imageData);
	}
}



class Point3
{
public:
	float x, y, z;
	void set(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}
	void set(Point3& p)
	{
		x = p.x; y = p.y; z = p.z;
	}
	Point3() { x = y = z = 0; }
	Point3(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}

};
class Color3
{
public:
	float r, g, b;
	void set(float red, float green, float blue)
	{
		r = red; g = green; b = blue;
	}
	void set(Color3& c)
	{
		r = c.r; g = c.g; b = c.b;
	}
	Color3() { r = g = b = 0; }
	Color3(float red, float green, float blue)
	{
		r = red; g = green; b = blue;
	}

};
class Point2
{
public:
	Point2() { x = y = 0.0f; } // constructor 1
	Point2(float xx, float yy) { x = xx; y = yy; } // constructor 2
	void set(float xx, float yy) { x = xx; y = yy; }
	float getX() { return x; }
	float getY() { return y; }
	void draw()		{
		glBegin(GL_POINTS);
		glVertex2f((GLfloat)x, (GLfloat)y);
		glEnd();
	}
private:
	float 	x, y;
};
class IntRect
{
public:
	IntRect() { l = 0; r = 100; b = 0; t = 100; } // constructor
	IntRect(int left, int right, int bottom, int top)
	{
		l = left; r = right; b = bottom; t = top;
	}
	void set(int left, int right, int bottom, int top)
	{
		l = left; r = right; b = bottom; t = top;
	}
	void draw(){
		glRecti(l, b, r, t);
		glFlush();
	} // draw this rectangle using OpenGL
	int  getWidth(){ return (r - l); }
	int  getHeight() { return (t - b); }
private:
	int	l, r, b, t;
};
class RealRect
{
public:
	RealRect() { l = 0; r = 100; b = 0; t = 100; } // constructor
	RealRect(float left, float right, float bottom, float top)
	{
		l = left; r = right; b = bottom; t = top;
	}
	void set(float left, float right, float bottom, float top)
	{
		l = left; r = right; b = bottom; t = top;
	}
	float  getWidth(){ return (r - l); }
	float  getHeight() { return (t - b); }
	void RealRect::draw(){
		glRectf(l, b, r, t);
		glFlush();
	};// draw this rectangle using OpenGL
private:
	float	l, r, b, t;
};
class Vector3
{
public:
	float	x, y, z;
	void set(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}
	void set(Vector3& v)
	{
		x = v.x; y = v.y; z = v.z;
	}
	void flip()
	{
		x = -x; y = -y; z = -z;
	}
	void normalize();
	Vector3() { x = y = z = 0; }
	Vector3(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}
	Vector3(Vector3& v)
	{
		x = v.x; y = v.y; z = v.z;
	}
	Vector3 cross(Vector3 b);
	float dot(Vector3 b);
};
class VertexID
{
public:
	int		vertIndex;
	int		colorIndex;
};
class Face
{
public:
	int		nVerts;
	VertexID*	vert;
	Vector3		facenorm;
	Face()
	{
		nVerts = 0;
		vert = NULL;
	}
	~Face()
	{
		if (vert != NULL)
		{
			delete[] vert;
			vert = NULL;
		}
		nVerts = 0;
	}
};
class Mesh
{
public:
	int		numVerts;
	Point3*		pt;

	int		numFaces;
	Face*		face;
public:
	Mesh()
	{
		numVerts = 0;
		pt = NULL;
		numFaces = 0;
		face = NULL;
		rotateX = 0;
		rotateY = 0;
		rotateZ = 0;
	}
	~Mesh()
	{
		if (pt != NULL)
		{
			delete[] pt;
		}
		if (face != NULL)
		{
			delete[] face;
		}
		numVerts = 0;
		numFaces = 0;
	}
	void DrawWireframe();
	void DrawColor();
	void Draw();
	Vector3 CalculateFacesNorm1(Point3 p1, Point3 p2, Point3 p3);
	void CalculateFacesNorm();
	void CreateCuboid(float	fSizeX, float fSizeY, float	fSizeZ);
	void CreateCylinder(int nSegment, float fHeight, float fRadius);
	void CreateDoubleCircle(float fR1, float fR2, float fLength, float fHeight, int nVertex);
	void CreateDiamond();
	float slideX, slideY, slideZ;
	float rotateX, rotateY, rotateZ;
	float scaleX, scaleY, scaleZ;
	void SetColor(int colorIdx);

};
class Camera
{
private:
	Point3      eye;
	Point3		look;
	Vector3      u, v, n;
	double       viewAngle, aspect, nearDist, farDist;
	void         setModelViewMatrix();
public:
	Camera(){};
	void	set(Point3 Eye, Point3 Look, Vector3 up);
	void	roll(float angle);
	void	pitch(float angle);
	void	yaw(float angle);
	void	slide(float delU, float delV, float delN);
	void	distanceY(float disY);
	void	distanceXZ(float height);
	void	moveUp(float step);
	void	moveDown(float step);
	void	moveNearFar(float step);
	void	moveRL(float step);
	void	setShape(float vAng, float asp, float nearD, float farD);
	float	radius(Point3 E);
};

Vector3 Vector3::cross(Vector3 b)
{
	Vector3 c(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
	return c;
}
float Vector3::dot(Vector3 b)
{
	return x*b.x + y*b.y + z*b.z;
}
void Vector3::normalize()
{
	float temp = sqrt(x*x + y*y + z*z);
	x = x / temp;
	y = y / temp;
	z = z / temp;
}

void Mesh::CalculateFacesNorm(){
	Vector3 v1,v2;
	Point3 p1, p2, p3;
	for (int i = 0; i < numFaces; i++){
		p1 = pt[face[i].vert[0].vertIndex];
		p2 = pt[face[i].vert[1].vertIndex];
		p3 = pt[face[i].vert[2].vertIndex];
		v1 = Vector3(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
		v2 = Vector3(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);
		face[i].facenorm = Vector3(v1.cross(v2));
		face[i].facenorm.normalize();
	}
}

void Mesh::CreateCuboid(float fSizeX, float fSizeY, float fSizeZ)
{
	int i;

	numVerts = 8;
	pt = new Point3[numVerts];
	pt[0].set(-fSizeX, fSizeY, fSizeZ);
	pt[1].set(fSizeX, fSizeY, fSizeZ);
	pt[2].set(fSizeX, fSizeY, -fSizeZ);
	pt[3].set(-fSizeX, fSizeY, -fSizeZ);
	pt[4].set(-fSizeX, -fSizeY, fSizeZ);
	pt[5].set(fSizeX, -fSizeY, fSizeZ);
	pt[6].set(fSizeX, -fSizeY, -fSizeZ);
	pt[7].set(-fSizeX, -fSizeY, -fSizeZ);


	numFaces = 6;
	face = new Face[numFaces];

	//Left face
	face[0].nVerts = 4;
	face[0].vert = new VertexID[face[0].nVerts];
	face[0].vert[0].vertIndex = 1;
	face[0].vert[1].vertIndex = 5;
	face[0].vert[2].vertIndex = 6;
	face[0].vert[3].vertIndex = 2;
	for (i = 0; i<face[0].nVerts; i++)
		face[0].vert[i].colorIndex = 0;

	//Right face
	face[1].nVerts = 4;
	face[1].vert = new VertexID[face[1].nVerts];
	face[1].vert[0].vertIndex = 0;
	face[1].vert[1].vertIndex = 3;
	face[1].vert[2].vertIndex = 7;
	face[1].vert[3].vertIndex = 4;
	for (i = 0; i<face[1].nVerts; i++)
		face[1].vert[i].colorIndex = 1;

	//top face
	face[2].nVerts = 4;
	face[2].vert = new VertexID[face[2].nVerts];
	face[2].vert[0].vertIndex = 0;
	face[2].vert[1].vertIndex = 1;
	face[2].vert[2].vertIndex = 2;
	face[2].vert[3].vertIndex = 3;
	for (i = 0; i<face[2].nVerts; i++)
		face[2].vert[i].colorIndex = 2;

	//bottom face
	face[3].nVerts = 4;
	face[3].vert = new VertexID[face[3].nVerts];
	face[3].vert[0].vertIndex = 7;
	face[3].vert[1].vertIndex = 6;
	face[3].vert[2].vertIndex = 5;
	face[3].vert[3].vertIndex = 4;
	for (i = 0; i<face[3].nVerts; i++)
		face[3].vert[i].colorIndex = 3;

	//near face
	face[4].nVerts = 4;
	face[4].vert = new VertexID[face[4].nVerts];
	face[4].vert[0].vertIndex = 4;
	face[4].vert[1].vertIndex = 5;
	face[4].vert[2].vertIndex = 1;
	face[4].vert[3].vertIndex = 0;
	for (i = 0; i<face[4].nVerts; i++)
		face[4].vert[i].colorIndex = 4;

	//Far face
	face[5].nVerts = 4;
	face[5].vert = new VertexID[face[5].nVerts];
	face[5].vert[0].vertIndex = 3;
	face[5].vert[1].vertIndex = 2;
	face[5].vert[2].vertIndex = 6;
	face[5].vert[3].vertIndex = 7;
	for (i = 0; i<face[5].nVerts; i++)
		face[5].vert[i].colorIndex = 5;

}
void Mesh::CreateCylinder(int nSegment, float fHeight, float fRadius){

	int i;
	int j;

	numVerts = nSegment + 1;
	pt = new Point3[numVerts * 2];
	//pt[0].set(fRadius*sin(2*PI/nSegment), fRadius*cos(2 * PI / nSegment),0);

	for (i = 0; i < nSegment; i++)
		pt[i].set(fRadius*sin(2 * i * PI / nSegment), 0, fRadius*cos(2 * i * PI / nSegment));

	for (i = nSegment; i < 2 * nSegment; i++)
		pt[i].set(fRadius*sin(2 * i * PI / nSegment), fHeight, fRadius*cos(2 * i * PI / nSegment));

	pt[2 * nSegment].set(0, 0, 0);
	pt[2 * nSegment + 1].set(0, fHeight, 0);

	numFaces = nSegment * 3;
	face = new Face[numFaces];

	//ve da giac cuoi cung
	/**/
	face[nSegment - 1].nVerts = 3;
	face[nSegment - 1].vert = new VertexID[face[nSegment - 1].nVerts];
	face[nSegment - 1].vert[0].vertIndex = nSegment - 1;
	face[nSegment - 1].vert[1].vertIndex = 0;
	face[nSegment - 1].vert[2].vertIndex = 2 * nSegment;

	face[2 * nSegment - 1].nVerts = 3;
	face[2 * nSegment - 1].vert = new VertexID[face[nSegment - 1].nVerts];
	face[2 * nSegment - 1].vert[0].vertIndex = 2 * nSegment - 1;
	face[2 * nSegment - 1].vert[1].vertIndex = nSegment;
	face[2 * nSegment - 1].vert[2].vertIndex = 2 * nSegment + 1;

	face[3 * nSegment - 1].nVerts = 4;
	face[3 * nSegment - 1].vert = new VertexID[face[3 * nSegment - 1].nVerts];
	face[3 * nSegment - 1].vert[0].vertIndex = nSegment - 1;
	face[3 * nSegment - 1].vert[1].vertIndex = 0;
	face[3 * nSegment - 1].vert[2].vertIndex = nSegment;
	face[3 * nSegment - 1].vert[3].vertIndex = 2 * nSegment - 1;


	for (i = 0; i < nSegment - 1; i++)
	{
		//mat tron 1
		face[i].nVerts = 3;
		face[i].vert = new VertexID[face[i].nVerts];

		face[i].vert[0].vertIndex = i;
		face[i].vert[1].vertIndex = i + 1;
		face[i].vert[2].vertIndex = 2 * nSegment;// tam duong tron 1

		for (j = 0; j < face[i].nVerts; j++)
			face[i].vert[j].colorIndex = i;

		//mat tron 2
		face[i + nSegment].nVerts = 3;
		face[i + nSegment].vert = new VertexID[face[i + nSegment].nVerts];

		face[i + nSegment].vert[0].vertIndex = i + nSegment;
		face[i + nSegment].vert[1].vertIndex = i + nSegment + 1;
		face[i + nSegment].vert[2].vertIndex = 2 * nSegment + 1;// tam duong tron 2

		for (j = 0; j < face[i + nSegment].nVerts; j++)
			face[i + nSegment].vert[j].colorIndex = i;

		//mat ben

		face[i + 2 * nSegment].nVerts = 4;
		face[i + 2 * nSegment].vert = new VertexID[face[i + 2 * nSegment].nVerts];

		face[i + 2 * nSegment].vert[0].vertIndex = i;
		face[i + 2 * nSegment].vert[1].vertIndex = i + 1;
		face[i + 2 * nSegment].vert[2].vertIndex = i + nSegment + 1;
		face[i + 2 * nSegment].vert[3].vertIndex = i + nSegment;

		for (j = 0; j < face[i + 2 * nSegment].nVerts; j++)
			face[i + 2 * nSegment].vert[j].colorIndex = i;

	}
}
void Mesh::CreateDoubleCircle(float fR1, float fR2, float fLength, float fHeight, int nVertex){
	int i;
	int j;

	numVerts = nVertex * 4;
	pt = new Point3[numVerts];

	for (i = 0; i < nVertex; i++) {
		// mat duoi
		pt[i].set(fR1*sin(i * PI / (nVertex - 1)), 0, fR1*cos(i * PI / (nVertex - 1)));
		pt[i + nVertex].set(fLength + fR2*sin(PI + i * PI / (nVertex - 1)), 0, fR2*cos(PI + i * PI / (nVertex - 1)));
		//mat tren
		pt[i + 2 * nVertex].set(fR1*sin(i * PI / (nVertex - 1)), 0.2, fR1*cos(i * PI / (nVertex - 1)));
		pt[i + 3 * nVertex].set(fLength + fR2*sin(PI + i * PI / (nVertex - 1)), 0.2, fR2*cos(PI + i * PI / (nVertex - 1)));
	}

	numFaces = 2 * nVertex+2;
	face = new Face[numFaces];

	// mat cuoi cua mat ben
	face[2 * nVertex - 1].nVerts = 4;
	face[2 * nVertex - 1].vert = new VertexID[face[2 * nVertex - 1].nVerts];

	face[2 * nVertex - 1].vert[0].vertIndex = 0;
	face[2 * nVertex - 1].vert[1].vertIndex = 2 * nVertex;
	face[2 * nVertex - 1].vert[2].vertIndex = 4 * nVertex - 1;
	face[2 * nVertex - 1].vert[3].vertIndex = 2 * nVertex - 1;

	for (j = 0; j < face[2 * nVertex - 1].nVerts; j++)
		face[2 * nVertex - 1].vert[j].colorIndex = 1;

	for (i = 0; i < 2 * nVertex - 1; i++)
	{
		//mat ben
		face[i].nVerts = 4;
		face[i].vert = new VertexID[face[i].nVerts];

		face[i].vert[0].vertIndex = i;
		face[i].vert[1].vertIndex = i + 1;
		face[i].vert[2].vertIndex = i + 2 * nVertex + 1;
		face[i].vert[3].vertIndex = i + 2 * nVertex;

		for (j = 0; j < face[i].nVerts; j++)
			face[i].vert[j].colorIndex = i;
	}

	face[2 * nVertex].nVerts = 2 * nVertex;
	face[2 * nVertex].vert = new VertexID[face[2 * nVertex].nVerts];

	//mat tren
	for (i = 0; i < 2 * nVertex; i++)
	{
		//face[2 * nVertex].vert[i].vertIndex = 4 * nVertex-1-i;
		face[2 * nVertex].vert[i].vertIndex = 2 * nVertex +i;
	}

	face[2 * nVertex+1].nVerts = 2 * nVertex;
	face[2 * nVertex+1].vert = new VertexID[face[2 * nVertex+1].nVerts];

	//mat duoi
	for (i = 0; i < 2 * nVertex; i++)
	{
		face[2 * nVertex+1].vert[i].vertIndex = 2*nVertex-1- i;
	}

	for (j = 0; j < face[2 * nVertex].nVerts; j++)
		face[2 * nVertex+1].vert[j].colorIndex = 2 * nVertex;
}
void Mesh::CreateDiamond(){
	int i;
	float x = sqrt(3.0) / 2;
	float fSize = 1;
	numVerts = 4;

	pt = new Point3[numVerts];
	pt[0].set(0, 0, 0);
	pt[1].set(-x, 0, 0.5);
	pt[2].set(0, 0, 1);
	pt[3].set(x, 0, 0.5);

	numFaces = 1;
	face = new Face[numFaces];

	face[0].nVerts = 4;
	face[0].vert = new VertexID[face[0].nVerts];
	face[0].vert[0].vertIndex = 0;
	face[0].vert[1].vertIndex = 1;
	face[0].vert[2].vertIndex = 2;
	face[0].vert[3].vertIndex = 3;

}
void Mesh::DrawWireframe()
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int f = 0; f < numFaces; f++)
	{
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++)
		{
			int		iv = face[f].vert[v].vertIndex;

			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

void Mesh::Draw() {
	for (int f = 0; f < numFaces; f++){
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++){
			int		iv = face[f].vert[v].vertIndex;
			glNormal3f(face[f].facenorm.x, face[f].facenorm.y, face[f].facenorm.z);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}


void Mesh::DrawColor()
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (int f = 0; f < numFaces; f++)
	{
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++)
		{
			int		iv = face[f].vert[v].vertIndex;
			int		ic = face[f].vert[v].colorIndex;

			//ic = f % COLORNUM;

			glColor3f(ColorArr[ic][0], ColorArr[ic][1], ColorArr[ic][2]);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}
void Mesh::SetColor(int colorIdx){
	for (int f = 0; f < numFaces; f++)
	{
		for (int v = 0; v < face[f].nVerts; v++)
		{
			face[f].vert[v].colorIndex = colorIdx;
		}
	}
}
void Camera::setModelViewMatrix()
{
	float         m[16];
	Vector3   eVec(eye.x, eye.y, eye.z);
	m[0] = u.x; m[4] = u.y; m[8] = u.z; m[12] = -eVec.dot(u);
	m[1] = v.x; m[5] = v.y; m[9] = v.z; m[13] = -eVec.dot(v);
	m[2] = n.x; m[6] = n.y; m[10] = n.z; m[14] = -eVec.dot(n);
	m[3] = 0; m[7] = 0; m[11] = 0; m[15] = 1.0;
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(m);
}
void Camera::set(Point3 Eye, Point3 Look, Vector3 up)
{
	eye.set(Eye);
	n.set(Eye.x - Look.x, eye.y - Look.y, eye.z - Look.z);
	u.set(up.cross(n));
	n.normalize();
	u.normalize();
	v.set(n.cross(u));
	setModelViewMatrix();
}
void Camera::slide(float delU, float delV, float delN)
{
	eye.x += delU*u.x + delV*v.x + delN*n.x;
	eye.y += delU*u.y + delV*v.y + delN*n.y;
	eye.z += delU*u.z + delV*v.z + delN*n.z;
	setModelViewMatrix();
}
void Camera::moveUp(float step){

	eye.y += step;
	set(eye, look, Vector3(0, 1, 0));
	setModelViewMatrix();
}
void Camera::moveDown(float step){

	eye.y -= step;
	set(eye, look, Vector3(0, 1, 0));
	setModelViewMatrix();
}
void Camera::moveNearFar(float step){
	eye.x += step*n.x;
	eye.z += step*n.z;
	set(eye, look, Vector3(0, 1, 0));
	setModelViewMatrix();
}
void Camera::moveRL(float step){
	angle += step;
	if (angle >= 360)
		angle = 360 - angle;
	float r = sqrt(eye.x*eye.x + eye.z*eye.z);
	eye.x = r*sin(angle*PI / 180);
	eye.z = r*cos(angle*PI / 180);
	set(eye, look, Vector3(0, 1, 0));
	setModelViewMatrix();
}
void Camera::roll(float angle)
{
	float cs = cos(3.14159265 / 180 * angle);
	float sn = sin(3.14159265 / 180 * angle);
	Vector3	t = u;
	u.set(cs*t.x - sn*v.x, cs*t.y - sn*v.y, cs*t.z - sn*v.z);
	v.set(sn*t.x + cs*v.x, sn*t.y + cs*v.y, sn*t.z + cs*v.z);
	setModelViewMatrix();
}
void Camera::pitch(float angle)
{
	float cs = cos(3.14159265 / 180 * angle);
	float sn = sin(3.14159265 / 180 * angle);
	Vector3	t = n;
	n.set(cs*t.x - sn*v.x, cs*t.y - sn*v.y, cs*t.z - sn*v.z);
	v.set(sn*t.x + cs*v.x, sn*t.y + cs*v.y, sn*t.z + cs*v.z);
	setModelViewMatrix();
}
void Camera::yaw(float angle)
{
	float cs = cos(3.14159265 / 180 * angle);
	float sn = sin(3.14159265 / 180 * angle);
	Vector3	t = u;
	u.set(cs*t.x - sn*n.x, cs*t.y - sn*n.y, cs*t.z - sn*n.z);
	n.set(sn*t.x + cs*n.x, sn*t.y + cs*n.y, sn*t.z + cs*n.z);
	setModelViewMatrix();
}
void Camera::setShape(float vAng, float asp, float nearD, float farD)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	viewAngle = vAng;
	aspect = asp;
	nearDist = nearD;
	farDist = farD;
	gluPerspective(viewAngle, aspect, nearDist, farDist);
	setModelViewMatrix();
}


int		screenWidth = 600;
int		screenHeight = 600;

bool	bWireFrame = false;

float	baseRadius = 0.6;
float	baseHeight = 0.2;
float	baseRotateStep = 5;

float	cylRadius = 0.3;
float	cylHeight = 1.0;
float   cylMaxScaleY = 2.0;
float	cylScaleStep = 0.05;

float	body1SizeX = 3;
float	body1SizeY = 0.2;
float	body1SizeZ = 0.8;

float	volumeStep = 5.0;

bool editMode = true;
bool turnOnOff = false;
bool fourCamera = false;
int count_auto = 0;


float camera_angle = 5, camera_angle_height = 10.0, camera_height = 0.1, camera_dis = 0.1;
int bodyColor;

Point3 eye_cam(0, 1, 10);

Mesh	base;
Mesh	cyl;
Mesh	cyl1;
Mesh	cyl2;
Mesh	cyl4;
Mesh	body1;
Mesh	body2;
Mesh	body3;
Mesh	body4;
Mesh	scroll;
Mesh	arm1;
Mesh	arm2;
Mesh    flor1;
Mesh    flor2;
Mesh    flor3;
Mesh	cuboid;
Camera	cam;

float d2r(int deg){
	return deg*PI / 180;
}

void drawAxis()
{
	glPushMatrix();

	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(-4, 0, 0);//x
	glVertex3f(4, 0, 0);

	glColor3f(0, 1, 0);
	glVertex3f(0, 0, 0);//y
	glVertex3f(0, 4, 0);

	glColor3f(0, 0, 1);
	glVertex3f(0, 0, -4);//z
	glVertex3f(0, 0, 4);
	glEnd();

	glPopMatrix();
}

void processTimer(int value){
	if (!editMode){
		float	fRInc;
		float	fAngle;
		float alpha, beta, gama;
		float a, b, c;
		base.rotateY += baseRotateStep / 5;
		if (base.rotateY > 360)
			base.rotateY -= 360;
		arm1.rotateY += baseRotateStep;
		if (arm1.rotateY > 360)
			arm1.rotateY -= 360;

		a = 3.0; b = 1.6;
		if (arm1.rotateY<180 & arm1.rotateY >= 0)
			alpha = arm1.rotateY;
		else
			alpha = 360 - arm1.rotateY;

		beta = asin(sin(d2r(alpha)) * b / a);
		gama = PI - d2r(alpha) - beta;
		c = sin(gama)*a / sin(d2r(alpha));

		if (arm1.rotateY >= 0 && arm1.rotateY <= 180)
			arm2.rotateY = -beta / PI * 180; // arm2
		else
			arm2.rotateY = beta / PI * 180; // arm2

		if (alpha != 0 && alpha != 180)
			scroll.slideX = 4.6 - c; // 4.65

		if (count_auto < 15){
			count_auto++;
			cyl.scaleY += baseHeight / 2;
			if (count_auto == 15){
				count_auto = 30;
			}
		}
		if (count_auto > 15) {
			count_auto--;
			cyl.scaleY -= baseHeight / 2;
			if (count_auto == 15){
				count_auto = 0;
			}
		}

		glutTimerFunc(50, processTimer, value);
		glutPostRedisplay();
	}
}

GLfloat ambientm[] = { 0.1, 0.1, 0.1, 1.0f };
GLfloat specularm[] = { 1.0, 1.0, 1.0, 1.0f };
GLfloat	shinessm = 100.8;

void setMaterial(float r, float g, float b){
	GLfloat diffusem[] = { r, g, b, 1.0f };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambientm);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffusem);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularm);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shinessm);
}

void drawBase()
{
	glPushMatrix();

	glTranslated(0, base.slideY, 0);
	glRotatef(base.rotateY, 0, 1, 0);

	if (bWireFrame)
		base.DrawWireframe();
	else{
		setMaterial(1.0, 0, 0);
		base.Draw();
	}

	glPopMatrix();
}

void drawCyl()
{
	glPushMatrix();
	glTranslated(0, cyl.slideY - baseHeight, 0);
	glScalef(1, cyl.scaleY + baseHeight * 5, 1);
	glRotatef(base.rotateY, 0, 1, 0);
	if (bWireFrame)
		cyl.DrawWireframe();
	else{
		setMaterial(0, 0, 1.0);
		cyl.Draw();
	}
	glPopMatrix();
}

void drawBody1()
{
	bodyColor = 4;
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(0, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 5, 0);
	if (bWireFrame)
		body1.DrawWireframe();
	else{
		setMaterial(0, 1.0, 0);
		body1.Draw();
	}
	glPopMatrix();
}

void drawBody2()
{
	bodyColor = 5;
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(3 * body1SizeX / 4, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 8, 0);
	//glScalef(1.0 / 4, 2, 1);
	if (bWireFrame)
		body2.DrawWireframe();
	else{
		setMaterial(0, 1.0, 0);
		body2.Draw();
	}

	glPopMatrix();
}

void drawBody3()
{
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(-1.0 * body1SizeX / 4, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 7, 2.0*body1SizeZ / 3);
	//glScalef(3.0 / 4, 1, 1.0 / 3);
	if (bWireFrame)
		body3.DrawWireframe();
	else{
		setMaterial(0, 1.0, 0);
		body3.Draw();
	}
	glPopMatrix();
}

void drawBody4()
{
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(-1.0 * body1SizeX / 4, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 7, -2.0*body1SizeZ / 3);
	//glScalef(3.0 / 4, 1, 1.0 / 3);

	if (bWireFrame)
		body4.DrawWireframe();
	else{
		setMaterial(0, 1.0, 0);
		body4.Draw();
	}
	glPopMatrix();
}

void drawWireArm1()
{
	glPushMatrix();
	glColor3f(1, 0, 0);

	glBegin(GL_LINE_STRIP);
	glVertex3f(0.4, 0, 0.4);
	glVertex3f(-1.8, 0, 0.4);
	glVertex3f(-1.8, 0, -0.4);
	glVertex3f(0.4, 0, -0.4);
	glVertex3f(0.4, 0, 0.4);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(0.4, 0.2, 0.4);
	glVertex3f(-1.8, 0.2, 0.4);
	glVertex3f(-1.8, 0.2, -0.4);
	glVertex3f(0.4, 0.2, -0.4);
	glVertex3f(0.4, 0.2, 0.4);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(0.4, 0, 0.4);
	glVertex3f(0.4, 0.2, 0.4);

	glVertex3f(-1.8, 0, 0.4);
	glVertex3f(-1.8, 0.2, 0.4);

	glVertex3f(-1.8, 0, -0.4);
	glVertex3f(-1.8, 0.2, -0.4);

	glVertex3f(0.4, 0, -0.4);
	glVertex3f(0.4, 0.2, -0.4);

	glEnd();
	
	glPopMatrix();
}

void drawWireArm2()
{
	glPushMatrix();
	glColor3f(1, 0, 0);

	glBegin(GL_LINE_STRIP);
	glVertex3f(0.2, 0, 0.2);
	glVertex3f(-3.2, 0, 0.2);
	glVertex3f(-3.2, 0, -0.2);
	glVertex3f(0.2, 0, -0.2);
	glVertex3f(0.2, 0, 0.2);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(0.2, 0.2, 0.2);
	glVertex3f(-3.2, 0.2, 0.2);
	glVertex3f(-3.2, 0.2, -0.2);
	glVertex3f(0.2, 0.2, -0.2);
	glVertex3f(0.2, 0.2, 0.2);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(0.2, 0, 0.2);
	glVertex3f(0.2, 0.2, 0.2);

	glVertex3f(-3.2, 0, 0.2);
	glVertex3f(-3.2, 0.2, 0.2);

	glVertex3f(-3.2, 0, -0.2);
	glVertex3f(-3.2, 0.2, -0.2);

	glVertex3f(0.2, 0, -0.2);
	glVertex3f(0.2, 0.2, -0.2);
	glEnd();

	glPopMatrix();
}

void drawWireScroll()
{
	glPushMatrix();
	glColor3f(1, 0, 0);

	glBegin(GL_LINE_STRIP);
	glVertex3f(0.3, -0.5, 0.3);
	glVertex3f(-0.3, -0.5, 0.3);
	glVertex3f(-0.3, -0.5, -0.3);
	glVertex3f(0.3, -0.5, -0.3);
	glVertex3f(0.3, -0.5, 0.3);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(0.3, 0.5, 0.3);
	glVertex3f(-0.3, 0.5, 0.3);
	glVertex3f(-0.3, 0.5, -0.3);
	glVertex3f(0.3, 0.5, -0.3);
	glVertex3f(0.3, 0.5, 0.3);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(0.3, -0.5, 0.3);
	glVertex3f(0.3, 0.5, 0.3);

	glVertex3f(-0.3, -0.5, 0.3);
	glVertex3f(-0.3, 0.5, 0.3);

	glVertex3f(-0.3, -0.5, -0.3);
	glVertex3f(-0.3, 0.5, -0.3);

	glVertex3f(0.3, -0.5, -0.3);
	glVertex3f(0.3, 0.5, -0.3);
	glEnd();


	glPopMatrix();
}

void drawArm1()
{

	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 1, -1);
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(3 * body1SizeX / 4, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 10, 0);
	glRotatef(arm1.rotateY, 0, 1, 0);
	//glScalef(1, 0.2, 1);
	if (bWireFrame){
		arm1.DrawWireframe();
	}
	else{
		if (check_arm1) {
			if (!check_shadow1 && !check_shadow2){
				drawWireArm1();
			}
			if (check_r) {
				arm1_r = 1.0; arm1_g = 0.0; arm1_b = 0.0;
			}
			else if (check_g) {
				arm1_r = 0.0; arm1_g = 1.0; arm1_b = 0.0;
			}
			else if (check_b) {
				arm1_r = 0.0; arm1_g = 0.0; arm1_b = 1.0;
			}

		}
		
		setMaterial(arm1_r, arm1_g, arm1_b);	
		//drawAxis();
		arm1.Draw();
	}
	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}



void drawArm2()
{
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 2, -1);
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(-2.4 + scroll.slideX, 0, 0);
	glRotatef(arm2.rotateY, 0, 1, 0);
	glTranslated(3, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 11, 0);
	//glScalef(1, 0.2, 1);
	if (bWireFrame) {
		arm2.DrawWireframe();
		
	}
	else{
		if (check_arm2) {
			if (!check_shadow1 && !check_shadow2){
				drawWireArm2();
			}
			if (check_r) {
				arm2_r = 1.0; arm2_g = 0.0; arm2_b = 0.0;
			}
			else if (check_g) {
				arm2_r = 0.0; arm2_g = 1.0; arm2_b = 0.0;
			}
			else if (check_b) {
				arm2_r = 0.0; arm2_g = 0.0; arm2_b = 1.0;
			}
		}
		setMaterial(arm2_r, arm2_g, arm2_b);
		
		arm2.Draw();
	}
	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawScroll()
{
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 3, -1);
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(-2.4 + scroll.slideX, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 8.5, 0);
	//glScalef(0.8, baseHeight*5, 0.8);
	if (bWireFrame)
		scroll.DrawWireframe();
	else{
		if (check_scroll) {
			if (!check_shadow1 && !check_shadow2){
				drawWireScroll();
			}
			if (check_r) {
				scroll_r = 1.0; scroll_g = 0.0; scroll_b = 0.0;
			}
			else if (check_g) {
				scroll_r = 0.0; scroll_g = 1.0; scroll_b = 0.0;
			}
			else if (check_b) {
				scroll_r = 0.0; scroll_g = 0.0; scroll_b = 1.0;
			}
		}
		
		setMaterial(scroll_r, scroll_g, scroll_b);
		scroll.Draw();
	}
	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawCyl1() // chot arm 1
{
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(3 * body1SizeX / 4, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 10, 0);
	//glScalef(1, baseHeight + 0.01, 1);
	if (bWireFrame)
		cyl1.DrawWireframe();
	else{
		setMaterial(1.0,0,0);
		cyl1.Draw();
	}
	glPopMatrix();
}

void drawCyl2() // chot arm2 trong
{
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(2.25, 0, 0);
	glRotatef(arm1.rotateY, 0, 1, 0);
	glTranslated(-1.6, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 10, 0);
	//glScalef(0.2, baseHeight * 2 + 0.01, 0.2);
	if (bWireFrame)
		cyl2.DrawWireframe();
	else{
		setMaterial(1.0, 0, 0);
		cyl2.Draw();
	}
	glPopMatrix();
}

void drawCyl4() // chot arm 2 ngoai
{
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(-2.4 + scroll.slideX, body1SizeY + cylHeight*cyl.scaleY + baseHeight * 10, 0);
	//glScalef(0.2, baseHeight * 2 + 0.01, 0.2);
	if (bWireFrame)
		cyl4.DrawWireframe();
	else{
		setMaterial(1.0, 0, 0);
		cyl4.Draw();
	}
	glPopMatrix();
}

void drawBrick(){
	glPushMatrix();
	for (int i = 0; i < 3; i++){
		glPushMatrix();
		glRotatef(120 * i, 0, 1, 0);
		switch (i)
		{
		case 0:
			if (bWireFrame)
				flor1.DrawWireframe();
			else
				flor1.Draw();
			break;
		case 1:
			if (bWireFrame)
				flor2.DrawWireframe();
			else
				flor2.Draw();
			break;
		case 2:
			if (bWireFrame)
				flor3.DrawWireframe();
			else
				flor3.Draw();
			break;
		default:
			break;
		}

		glPopMatrix();
	}
	glPopMatrix();
}

void drawFloor(){
	glDisable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);
	glColor4f(1.0, 1.0, 1.0, 0.40);  /* 40% dark red floor color */

	
	
	glBindTexture(GL_TEXTURE_2D, floorTex.texID);
	glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(-10.0f, 0.0f, -10.0f);
		glTexCoord2f(1.0, 0.0);
		glVertex3f(-10.0f, 0.0f, 10.0f);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(10.0f, 0.0f, 10.0f);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(10.0f, 0.0f, -10.0f);
	glEnd();
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
}

void shadowMatrix(GLfloat m[4][4], GLfloat plane[4], GLfloat light[4])
{
	GLfloat dot = plane[0] * light[0] + plane[1] * light[1] +
		plane[2] * light[2] + plane[3] * light[3];
	m[0][0] = dot - light[0] * plane[0];
	m[1][0] = -light[0] * plane[1];
	m[2][0] = -light[0] * plane[2];
	m[3][0] = -light[0] * plane[3];
	m[0][1] = -light[1] * plane[0];
	m[1][1] = dot - light[1] * plane[1];
	m[2][1] = -light[1] * plane[2];
	m[3][1] = -light[1] * plane[3];
	m[0][2] = -light[2] * plane[0];
	m[1][2] = -light[2] * plane[1];
	m[2][2] = dot - light[2] * plane[2];
	m[3][2] = -light[2] * plane[3];
	m[0][3] = -light[3] * plane[0];
	m[1][3] = -light[3] * plane[1];
	m[2][3] = -light[3] * plane[2];
	m[3][3] = dot - light[3] * plane[3];
}

void drawAll(){
	//drawAxis();
	drawBase();
	drawCyl();
	drawBody1();
	drawBody2();
	drawBody3();
	drawBody4();
	drawArm1();
	drawCyl1();
	drawArm2();
	drawScroll();// tru truot
	drawCyl2();
}

void myKeyboard(unsigned char key, int x, int y)
{
	float	fRInc;
	float	fAngle;
	float alpha, beta, gama;
	float a, b, c;
	switch (key)
	{
	case '1':
		base.rotateY += baseRotateStep;
		if (base.rotateY > 360)
			base.rotateY -= 360;
		break;
	case '2':
		base.rotateY -= baseRotateStep;
		if (base.rotateY < 0)
			base.rotateY += 360;
		break;
	case '3':
		arm1.rotateY += baseRotateStep;
		if (arm1.rotateY > 360)
			arm1.rotateY -= 360;

		a = 3.0; b = 1.6;
		if (arm1.rotateY<180 & arm1.rotateY >= 0)
			alpha = arm1.rotateY;
		else
			alpha = 360 - arm1.rotateY;

		beta = asin(sin(d2r(alpha)) * b / a);
		gama = PI - d2r(alpha) - beta;
		c = sin(gama)*a / sin(d2r(alpha));

		if (arm1.rotateY >= 0 && arm1.rotateY <= 180)
			arm2.rotateY = -beta / PI * 180; // arm2
		else
			arm2.rotateY = beta / PI * 180; // arm2

		if (alpha != 0 && alpha != 180)
			scroll.slideX = 4.6 - c; // 4.65

		break;
	case '4':
		arm1.rotateY -= baseRotateStep;
		if (arm1.rotateY > 360)
			arm1.rotateY -= 360;

		a = 3.0; b = 1.6;
		if (arm1.rotateY<180 & arm1.rotateY >= 0)
			alpha = arm1.rotateY;
		else
			alpha = 360 - arm1.rotateY;

		beta = asin(sin(d2r(alpha)) * b / a);
		gama = PI - d2r(alpha) - beta;
		c = sin(gama)*a / sin(d2r(alpha));

		if (arm1.rotateY >= 0 && arm1.rotateY <= 180)
			arm2.rotateY = -beta / PI * 180; // arm2
		else
			arm2.rotateY = beta / PI * 180; // arm2

		if (alpha != 0 && alpha != 180)
			scroll.slideX = 4.6 - c; // 4.65

		break;
	case 'l':
	case 'L':
		cyl.scaleY += baseHeight;
		break;
	case 'x':
	case 'X':
		cyl.scaleY -= baseHeight;
		break;
	case 'w':
	case 'W':
		bWireFrame = !bWireFrame;
		break;
	case '+':
		cam.moveNearFar(camera_dis);
		break;
	case '-':
		cam.moveNearFar(-camera_dis);
		break;
	case 'a':
	case 'A':
		editMode = !editMode;
		glutTimerFunc(50, processTimer, 1);
		break;

	case 'v':
	case 'V':
		fourCamera = !fourCamera;
		break;
	case 'd':
	case 'D':
		turnOnOff = !turnOnOff;
		if (turnOnOff){
			GLfloat	lightDiffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			GLfloat	lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			GLfloat	lightAmbient[] = { 0.4f, 0.4f, 0.4f, 1.0f };
			GLfloat light_position2[] = { -6.0f, 6.0f, 6.0f, 0.0f };
			GLfloat floor_position[4] = { 0.0, 1.0, 0.0, -0.1 };
			glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
			glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse);
			glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient);
			glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular);
			glEnable(GL_LIGHT1);

			

		}
		else{
			glDisable(GL_LIGHT1);
		}
		break;
	case 'r':
	case 'R':
		check_r = true;
		check_b = false;
		check_g = false;
		break;
	case 'b':
	case 'B':
		check_b = true;
		check_r = false;
		check_g = false;
		break;
	case 'g':
	case 'G':
		check_g = true;
		check_r = false;
		check_b = false;
		break;
	}
	glutPostRedisplay();
}



void myKeyboard2(int key, int x, int y)
{
	switch (key)
	{
		break;
	case GLUT_KEY_UP:
		cam.moveUp(camera_height);
		break;
	case GLUT_KEY_DOWN:
		cam.moveDown(camera_height);
		break;
	case GLUT_KEY_RIGHT:
		cam.moveRL(camera_angle);
		break;
	case GLUT_KEY_LEFT:
		cam.moveRL(-camera_angle);
		break;
	}

	glutPostRedisplay();
}
int window_width;
int window_height;

void clickMouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		window_width = glutGet(GLUT_WINDOW_WIDTH);
		window_height = glutGet(GLUT_WINDOW_HEIGHT);

		GLuint index;

		glReadPixels(x, window_height - y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &index);
		if (index == 1){
			check_arm1 = true;
			check_arm2 = false;
			check_scroll = false;
			check_r = false;
			check_g = false;
			check_b = false;
		}
		else if (index == 2){
			check_arm2 = true;
			check_arm1 = false;
			check_scroll = false;
			check_r = false;
			check_g = false;
			check_b = false;
		}
		else if (index = 3) {
			check_arm2 = false;
			check_arm1 = false;
			check_scroll = true;
			check_r = false;
			check_g = false;
			check_b = false;
		}
		else {
			check_arm1 = false;
			check_arm2 = false;
			check_scroll = false;
			check_r = false;
			check_g = false;
			check_b = false;
		}
		glutPostRedisplay();
	}

	
	
}

void myDisplay()
{


	GLfloat light_position1[4] = { 6.0f, 6.0f, 6.0f, 0.0f };
	GLfloat light_position2[] = { -6.0f, 6.0f, 6.0f, 0.0f };
	GLfloat floor_position[4] = { 0.0, 1.0, 0.0, -0.1 };

	glClearStencil(0); // this is the default value
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);	
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	if (fourCamera == false) {


		//shadow--------
		glPushMatrix();
		GLfloat matrix[4][4];
		shadowMatrix(&matrix[0], floor_position, light_position1);
		glMultMatrixf(&matrix[0][0]);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 2.0);
		glDisable(GL_LIGHTING);
		glColor3f(0.5, 0.5, 0.5);
		check_shadow1 = true;
		drawAll();
		check_shadow1 = false;
		glDisable(GL_POLYGON_OFFSET_FILL);
		glEnable(GL_LIGHTING);
		glPopMatrix();
		//--------------
		
		if (turnOnOff){
			glPushMatrix();
			GLfloat matrix2[4][4];
			shadowMatrix(&matrix2[0], floor_position, light_position2);
			glMultMatrixf(&matrix2[0][0]);
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0, 2.0);
			glDisable(GL_LIGHTING);
			glColor3f(0.5, 0.5, 0.5);
			check_shadow2 = true;
			drawAll();
			check_shadow2 = false;
			glDisable(GL_POLYGON_OFFSET_FILL);
			glEnable(GL_LIGHTING);
			glPopMatrix();
		}
		
		//phan chieu-----------
		glPushMatrix();
		glTranslatef(0, 0.2, 0);
		glScalef(1.0, -1.0, 1.0);
		drawAll();
		glPopMatrix();

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		drawFloor();
		glDisable(GL_BLEND);
		//--------------------
		
		drawAll();

		//drawWire(1, 1, 1);

		
		
		

	}
	else {
		for (int i = 0; i < 4; i++){

			switch (i)
			{
			case 0:
				cam.set(Point3(0, 10, 0.01), Point3(0, 0, 0), Vector3(0, 1, 0));
				glViewport(0, 0, 320, 240);

				break;
			case 3:
				cam.set(Point3(0, 1, 8), Point3(0, 0, 0), Vector3(0, 1, 0));
				glViewport(320, 0, 320, 240);
				break;
			case 2:

				cam.set(Point3(0, 0, 8), Point3(0, 0, 0), Vector3(0, 1, 0));
				glViewport(0, 240, 320, 240);
				break;
			case 1:
				cam.set(Point3(8, 0, 0), Point3(0, 0, 0), Vector3(0, 1, 0));
				glViewport(320, 240, 320, 240);

				break;
			default:
				break;
			}

			drawAxis();
			drawBase();
			drawCyl();

			drawBody1();
			drawBody2();
			drawBody3();
			drawBody4();

			drawArm1();
			drawCyl1();
			drawArm2();
			drawScroll();// tru truot
			drawCyl2();
			drawCyl4();
		}
	}

	glFlush();
	glutSwapBuffers();




}

void myInit()
{
	float	fHalfSize = 10000;
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
	loadTextures();

	GLfloat	lightDiffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat	lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat	lightAmbient[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat light_position1[] = { 6.0f, 6.0f, 6.0f, 0.0f };
	GLfloat floor_position[] = { 1.0, 0.0, 1.0, 0.0 };
	GLfloat light_position2[] = { -6.0f, 6.0f, 6.0f, 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glColor3f(0.0f, 0.0f, 0.0f);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-fHalfSize, fHalfSize, -fHalfSize, fHalfSize, -10000, 10000);
}

int main(int argc, char* argv[]){
	cout << "1, 2: Rotate the base" << endl;
	cout << "3, 4: Rotate the arm" << endl;
	cout << "L, l: Cylinder up" << endl;
	cout << "X, x: Cylinder down" << endl;
	cout << "W, w: Switch between wireframe and solid mode" << endl;
	cout << "V, v: to switch between 1 and 4 views." << endl;
	cout << "A, a: Turn on/of animation." << endl;
	cout << "+	: to increase camera distance." << endl;
	cout << "-	: to decrease camera distance." << endl;
	cout << "up arow  : to increase camera height." << endl;
	cout << "down arow: to decrease camera height." << endl;
	cout << "<-	      : to rotate camera clockwise." << endl;
	cout << "->	      : to rotate camera counterclockwise." << endl;




	glutInit(&argc, argv); //initialize the tool kit
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_STENCIL);//set the display mode
	glutInitWindowSize(640, 480); //set window size
	glutInitWindowPosition(500, 150); // set window position on screen
	glutCreateWindow("Assignment2 - Nguyen Thanh Hai(51301051)"); // open the screen window



	base.CreateCylinder(20, baseHeight, baseRadius);
	base.CalculateFacesNorm();
	base.SetColor(2);
	base.slideY = baseHeight / 2.0;

	cyl.CreateCylinder(20, cylHeight, cylRadius);
	cyl.CalculateFacesNorm();
	cyl.SetColor(0);
	cyl.slideY = cylHeight / 2.0;

	//glScalef(1, baseHeight + 0.01, 1);
	cyl1.CreateCylinder(20, baseHeight + 0.01, cylRadius);
	cyl1.CalculateFacesNorm();

	//glScalef(0.2, baseHeight * 2 + 0.01, 0.2);
	cyl2.CreateCylinder(20, baseHeight*2 + 0.01, cylRadius*0.2);
	cyl2.CalculateFacesNorm();

	//glScalef(0.2, baseHeight * 2 + 0.01, 0.2);
	cyl4.CreateCylinder(20, baseHeight*2 + 0.01, cylRadius*0.2);
	cyl4.CalculateFacesNorm();

	body1.CreateCuboid(body1SizeX, body1SizeY, body1SizeZ);
	body1.CalculateFacesNorm();
	//body1.SetColor(4);

	//glScalef(1.0 / 4, 2, 1);
	body2.CreateCuboid(body1SizeX/4, body1SizeY*2, body1SizeZ);
	body2.CalculateFacesNorm();
	//body2.SetColor(5);

	//glScalef(3.0 / 4, 1, 1.0 / 3);
	body3.CreateCuboid(3*body1SizeX/4, body1SizeY, body1SizeZ/3);
	body3.CalculateFacesNorm();
	//body3.SetColor(3);

	//glScalef(3.0 / 4, 1, 1.0 / 3);
	body4.CreateCuboid(3*body1SizeX/4, body1SizeY, body1SizeZ/3);
	body4.CalculateFacesNorm();
	//body4.SetColor(2);

	scroll.CreateCuboid(0.8 / 3, 0.5, 0.8 / 3);
	scroll.CalculateFacesNorm();
	//scroll.SetColor(0);
	//glScalef(0.8, baseHeight*5, 0.8);

	//glScalef(1, 0.2, 1);
	arm1.CreateDoubleCircle(body1SizeZ / 2, body1SizeZ / 4, -1.6, 1.0, 18);
	arm1.CalculateFacesNorm();
	//arm1.SetColor(6);

	//glScalef(1, 0.2, 1);
	arm2.CreateDoubleCircle(body1SizeZ / 4, body1SizeZ / 4, -3.0, 1.0, 18);
	arm2.CalculateFacesNorm();
	//arm2.SetColor(2);

	myInit();
	glutMouseFunc(clickMouse);
	glutKeyboardFunc(myKeyboard);
	glutSpecialFunc(myKeyboard2);
	
	
	glutDisplayFunc(myDisplay);
	glutTimerFunc(50, processTimer, 1);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glColor3f(0.0f, 0.0f, 0.0f);

	cam.set(eye_cam, Point3(0, 0, 0), Vector3(0, 1, 0));
	cam.setShape(40.0f, 64.0f / 48.0f, 0.5f, 500.0f);
	glutMainLoop();
	return 0;
}
