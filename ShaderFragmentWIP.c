#version 120

///////////////////////////   Type Definitions   //////////////////////////////////

varying vec3 fNormal;
varying float time0_1;
varying float time0_2PI;
uniform mat4 modelViewMatrix;
varying vec2 Position;


vec4 Black =vec4(0.0,0.0,0.0,0.0);
vec4 Green =vec4(0.0,1.0,0.0,1.0);
vec4 Blue =vec4(0.0,0.0,1.0,1.0);
vec4 Yellow =vec4(0.0,0.1,1.0,1.0);
vec4 Red =vec4(1.0,0.0,0.0,1.0);
vec4 White =vec4(1.0,1.0,1.0,1.0);
vec4 Grey =vec4(0.5,0.5,0.5,1.0);

vec3 X_Axis = vec3(1.0,0,0);
vec3 Y_Axis = vec3(0,1.0,0);
vec3 Z_Axis = vec3(0,0,1.0);


/////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////////   Defines  //////////////////////////////////////

#define INTERPOLATION_STEPS 1
#define MAX_RANGE_PARTICLES 15.0
//Constants
#define PI             3.14159
#define PI2             (PI*2.0)
#define INT_MAX          ((!0)-1)
#define INT_MIN          (!0)
#define min(x,y)       (x > y ? x : y)

//debugInfo
#define isParticleNearPosition(PARTICLE,POSITION,TOLERANCE)  (distance(PARTICLE + Center.Pos,POSITION) < TOLERANCE) 
//Helper functions
#define RAND(X,Y)        (fract(sin(dot(X + Y ,vec2(12.9898,78.233)))*43758.5453))
#define MATH_MAX(A,B)     = (A > B)? A: B;


#define norm(v)           sqrt(dot(v,v))     // norm = length of  vector
#define dist(u,v)         norm(u-v)          // distance = norm of difference
#define GRADIENT_SIZE    2
#define GRADIENTSIZE_FLOAT 2.0



#define access(MAT, ROW, COLUMN) MAT.a[((ROW)*8)+ (COLUMN)]

#define mat8RowMulRec(arr, row, rec)    arr[row*8+0]* rec.fLow.x +\
	arr[row*8+1]* rec.fUp.x   +\
	arr[row*8+2]* rec.rUp.x   +\
	arr[row*8+3]* rec.rLow.x   +\
	arr[row*8+4]* rec.fLow.y +\
	arr[row*8+5]* rec.fUp.y   +\
	arr[row*8+6]* rec.rUp.y   +\
	arr[row*8+7]* rec.rLow.y   

#define FRONT true
#define REAR false

///////////////////////////DataStructures  ///////////////////////////////////////
//Describes a normsquare Colour Gradient by a front and rear Gradient
struct colSampler{
	vec4 frontGradient[GRADIENT_SIZE ];
	vec4 rearGradient [GRADIENT_SIZE ];
};

//Works as a mask, allowing to precomp, wether a rotation is needed
struct rectangle{
	vec2 fLow;
	vec2 fUp;
	vec2 rLow;
	vec2 rUp;   
};

struct mat8 {
	float a[64];
};


///////////////////////////   Globals   //////////////////////////////////////

rectangle  warpedRectangle;
rectangle   normSquare;
colSampler  normSquareGradient;

///////////////////////////   Mathfunctions   //////////////////////////////////////
mat8 initializeProjectionMatrice(rectangle rec) {
	mat8 projectionMatrix;
	projectionMatrix.a[0] =    0;
	projectionMatrix.a[1] =    0;
	projectionMatrix.a[2] =    1;
	projectionMatrix.a[3] =    0;
	projectionMatrix.a[4] =    0;
	projectionMatrix.a[5] =    0;
	projectionMatrix.a[6] =    0;
	projectionMatrix.a[7] =    0;

	projectionMatrix.a[8+0] =    1;
	projectionMatrix.a[8+1] =    0;
	projectionMatrix.a[8+2] =    1;
	projectionMatrix.a[8+3] =    0;
	projectionMatrix.a[8+4] =    0;
	projectionMatrix.a[8+5] =    0;
	projectionMatrix.a[8+6] =    -1*rec.fUp.x;
	projectionMatrix.a[8+7] =    0;

	projectionMatrix.a[16+0] =      1;
	projectionMatrix.a[16+1] =      1;
	projectionMatrix.a[16+2] =      1;
	projectionMatrix.a[16+3] =      0;
	projectionMatrix.a[16+4] =      0;
	projectionMatrix.a[16+5] =      0;
	projectionMatrix.a[16+6] =      -1*rec.rUp.x;
	projectionMatrix.a[16+7] =      -1*rec.rUp.x;

	projectionMatrix.a[24+0] =     0;
	projectionMatrix.a[24+1] =     1;
	projectionMatrix.a[24+2] =     1;
	projectionMatrix.a[24+3] =     0;
	projectionMatrix.a[24+4] =     0;
	projectionMatrix.a[24+5] =     0;
	projectionMatrix.a[24+6] =     0;
	projectionMatrix.a[24+7] =     -1*rec.rLow.x;

	projectionMatrix.a[32+0] =     0;
	projectionMatrix.a[32+1] =     0;
	projectionMatrix.a[32+2] =     0;
	projectionMatrix.a[32+3] =     0;
	projectionMatrix.a[32+4] =     0;
	projectionMatrix.a[32+5] =     1;
	projectionMatrix.a[32+6] =     0;
	projectionMatrix.a[32+7] =     0; 

	projectionMatrix.a[40+0] =     0;
	projectionMatrix.a[40+1] =     0;
	projectionMatrix.a[40+2] =     0;
	projectionMatrix.a[40+3] =     1;
	projectionMatrix.a[40+4] =     0;
	projectionMatrix.a[40+5] =     1;
	projectionMatrix.a[40+6] =     -1*rec.fUp.y;
	projectionMatrix.a[40+7] =     0;

	projectionMatrix.a[48+0] =     0;
	projectionMatrix.a[48+1] =     0;
	projectionMatrix.a[48+2] =     0;
	projectionMatrix.a[48+3] =     1;
	projectionMatrix.a[48+4] =     1;
	projectionMatrix.a[48+5] =     1;
	projectionMatrix.a[48+6] =     -1*rec.rUp.y;
	projectionMatrix.a[48+7] =     -1*rec.rUp.y;

	projectionMatrix.a[56+0] =     0;
	projectionMatrix.a[56+1] =     0;
	projectionMatrix.a[56+2] =     0;
	projectionMatrix.a[56+3] =     0;
	projectionMatrix.a[56+4] =     1;
	projectionMatrix.a[56+5] =     1;
	projectionMatrix.a[56+6] =     0;
	projectionMatrix.a[56+7] =     -1*rec.rLow.y;


	return projectionMatrix;
}      


/*For calculating Determinant of the Matrix 
TODO: Make non recursive
*/

float determinantDecomposedMatrice(mat8 diagonalDeterminant) {
	return (diagonalDeterminant.a[0] *
	diagonalDeterminant.a[9] *
	diagonalDeterminant.a[18] *
	diagonalDeterminant.a[36] *
	diagonalDeterminant.a[45] *
	diagonalDeterminant.a[54] *
	diagonalDeterminant.a[63] );
}

//Based upon the crout algo
float determinant(mat8 A, int n) {
	int i, j, k;
	float sum = float( 0.0);
	mat8 L;
	mat8 U;

	for (i = 0; i < n; i++) {
		access(U,i,i) = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + access(L,i,k) * access(U,k,j);   
			}
			access(L,i,j) = access(A,i,j) - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + access(L,j,k) * access(U,k,i);
			}
			if (access(L,j,j) == 0) {
				access(L,j,j)= 0.00001;
			}
			access(U,j,i) = (access(A,j,i) - sum) / access(L,j,j);
		}
	}
	return determinantDecomposedMatrice(U) * determinantDecomposedMatrice(L);
}


mat8 transpose(mat8 num, mat8 fac, int r);

//Returns the Inverse
mat8 cofactor(mat8 num, int f)
{
	mat8 b;
	mat8 fac;

	int p, q, m, n, i, j;

	for (q = 0;q < f; q++)
	{
		for (p = 0;p < f; p++)
		{
			m = 0;
			n = 0;

			for (i = 0;i < f; i++)
			{
				for (j = 0;j < f; j++)
				{
					if (i != q && j != p)
					{
						access(b,m,n) = access(num,i,j);
						if (n < (f - 2)){ n++;}   else { n = 0; m++; }
					} 
				} 
			} 
			access(fac,q,p) = pow(-1, q + p) * determinant(b, f - 1);
		}
	}
	return transpose(num, fac, f);
}

/*Finding transpose of matrix*/ 
mat8 transpose(mat8 num, mat8 fac, int r)
{

	int i, j;
	mat8 b, inverse;
	float d;

	for (i = 0;i < r; i++)
	{
		for (j = 0;j < r; j++){
			access(b,i,j) = access(fac,j,i);
		}
	}
	d = determinant(num, r);

	for (i = 0;i < r; i++)
	{
		for (j = 0;j < r; j++)
		{
			access(inverse,i,j) = access(b,i,j) / d;
		}
	}
	return inverse;
}
///////////////////////////   Helperfunctions   //////////////////////////////////////

//retrieves a Colour Gradient from a array in a coloursampler -
vec4 getGradientAt(float fIndex, colSampler gradient, bool front) {
	float deNormalizedIndex = float(fIndex * GRADIENTSIZE_FLOAT);
	int   index = int(max(0,min(GRADIENT_SIZE,int(deNormalizedIndex))));
	int   upIndex = int(max(0,min(GRADIENT_SIZE,int(deNormalizedIndex)+1)));

	if (front == true) 
		return mix(gradient.frontGradient[index],gradient.frontGradient[upIndex], 
		fIndex);
	}
	else {
		return mix(gradient.rearGradient[index],gradient.rearGradient[upIndex], 
		fIndex);
	}
	return Black;
}


// determinates wether a point is within a triangle
bool PointInTriangle(vec2 a, vec2 b, vec2 c, vec2 s)
{
	float as_x = s.x-a.x;
	float as_y = s.y-a.y;

	bool s_ab = (b.x-a.x)*as_y-(b.y-a.y)*as_x > 0.0;

	if((c.x-a.x)*as_y-(c.y-a.y)*as_x > 0.0 == s_ab) return false;

	if((c.x-b.x)*(s.y-b.y)-(c.y-b.y)*(s.x-b.x) > 0.0 != s_ab) return false;

	return true;
}


// dist_Point_to_Segment(): get the distance of a point to a segment
//     Input:  a Point P and a Segment S (in any dimension)
//     Return: the shortest distance from P to S
float minimum_distance(vec2 P0, vec2 P1, vec2 P)
{
	vec2  v = P1 - P0;
	vec2  w = P - P0;

	float c1 = dot(w,v);
	if ( c1 <= 0 )
	return dist(P, P0);

	float c2 = dot(v,v);
	if ( c2 <= c1 )
	return dist(P, P1);


	vec2 Pb = P0 + (c1 / c2)*v;
	return dist(P, Pb);
}

// calculates the average of a rectangle edgepoints -
// does not  work for rectangles warped beyond a point
vec2 getRectangleCenter(rectangle rec){
	vec2 midValue = vec2(0,0);

	midValue +=rec.fUp;
	midValue +=rec.fLow;
	midValue +=rec.rUp;
	midValue +=rec.rLow;
	return midValue/4.0;
}


//debugPrintBorders of Quads
bool PointOnBorder(rectangle rec, vec2 pPos,float tol)
{
	vec2 offSet = getRectangleCenter(rec);
	rec.fUp +=offSet;
	rec.fLow +=offSet;
	rec.rUp +=offSet;
	rec.rLow +=offSet;
	if (minimum_distance(rec.fUp,rec.fLow,pPos) < tol) return true;
	if (minimum_distance(rec.fLow,rec.rLow,pPos) < tol) return true;
	if (minimum_distance(rec.rUp,rec.rLow,pPos) < tol) return true;
	if (minimum_distance(rec.rUp,rec.fUp,pPos) < tol) return true;
	return false;
}

//composes colours of a maped down point
vec4 recomposeColour(vec2 normSquareCoords, colSampler targetCol) {
	vec4 frontRowCol = getGradientAt(normSquareCoords.x,targetCol, FRONT);
	vec4 rearRowCol = getGradientAt(normSquareCoords.x,targetCol, REAR);

	return mix(frontRowCol, rearRowCol, normSquareCoords.y);      
}

// classic matrice * point multiplication                   
mat3 mulMatWithPoint(mat8 invertedProj, rectangle rec) {
	mat3 cof;  

	cof[0][0] = mat8RowMulRec(invertedProj.a, 0, rec);//a
	cof[0][1] = mat8RowMulRec(invertedProj.a, 1, rec);//b
	cof[0][2] = mat8RowMulRec(invertedProj.a, 2, rec);//c
	cof[1][0] = mat8RowMulRec(invertedProj.a, 3, rec);//d
	cof[1][1] = mat8RowMulRec(invertedProj.a, 4, rec);//e
	cof[1][2] = mat8RowMulRec(invertedProj.a, 5, rec);//f
	cof[2][0] = mat8RowMulRec(invertedProj.a, 6, rec);//g
	cof[2][1] = mat8RowMulRec(invertedProj.a, 7, rec);//h

	return cof;
}

//////////////////////////////////////////////////////////////////////////////////

vec4 drawRectangle(vec2 Position, rectangle rec, vec4 orgCol, vec4 targetCol) {
	if (PointInTriangle(rec.fLow,rec.fUp,rec.rUp, Position) ||
			PointInTriangle(rec.rUp,rec.rLow,rec.fLow, Position)) {
		return targetCol;
	}

	return orgCol;
}

bool isInRectangle(rectangle rec, vec2 position) {
	return  PointInTriangle(rec.fLow,rec.fUp,rec.rUp, position) ||
	PointInTriangle(rec.rUp,rec.rLow,rec.fLow, position);
}

//applys the derived projection coefficients
vec2 applyProjection(vec2 orgPos, mat3 cof) {
	vec2 result;
	result.x = (orgPos.x * cof[0][0] + orgPos.y * cof[0][1] + cof[0][2])/
	( orgPos.x *cof[2][0] + orgPos.y*cof[2][1] + 1);

	result.y = (orgPos.x * cof[1][0] + orgPos.y * cof[1][1] + cof[1][2])/
	( orgPos.x * cof[2][0] + orgPos.y* cof[2][1] + 1);

	return result; 
}

//computes the coefficients named in the paper
mat3 getCoefficients(rectangle rec) {

	mat8 org = initializeProjectionMatrice(rec);

	mat8 invertedProj = cofactor(org, 8);   

	return mulMatWithPoint(invertedProj, rec);
}

//Rectangle is in Square
vec4 drawWarpedRectangle(vec2 pos, rectangle rec, vec4 orgCol, colSampler targetCol) {
	if (isInRectangle(rec, pos)) {
		vec4 resultCol;
		
		// computate the coefficents of the projection matrice
		mat3 cof = getCoefficients(rec);
		
		// apply the projection matrice on the position
		vec2 normSquareCoords = applyProjection(pos, cof);
		
		// recompose color      
		resultCol = recomposeColour(normSquareCoords, targetCol);
		
		return resultCol;
	}

	return orgCol;
}

//initilizes testgeometry TODO THROW_AWAY
void init() {
	warpedRectangle.fLow = vec2(-3,-4);
	warpedRectangle.fUp = vec2(-4,5);
	warpedRectangle.rUp = vec2(6,7);
	warpedRectangle.rLow = vec2(6,-7);

	normSquare.fLow = vec2(-1,-1);
	normSquare.fUp = vec2(-1,1);
	normSquare.rUp = vec2(1,1);
	normSquare.rLow = vec2(1,-1);

	for (int i= 0;i < GRADIENT_SIZE; ++i) {
		normSquareGradient.frontGradient[i]= mix(Red,Blue,float(i)/GRADIENT_SIZE);
	}

	for (int i= 0;i < GRADIENT_SIZE; ++i) {
		normSquareGradient.rearGradient[i]= mix(Green,Yellow,float(i)/GRADIENT_SIZE);
	}
}

void main(void)
{      
	init();
	vec4 orgCol= Grey;
	orgCol = drawWarpedRectangle(Position, warpedRectangle, orgCol, normSquareGradient);
	orgCol = drawRectangle(Position, normSquare, orgCol, Red);

	gl_FragColor= orgCol;        

}


