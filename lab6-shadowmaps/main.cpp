
#include <windows.h>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#include "GL/glew.h"
#include "GL/glut.h"
#include "IL/ilut.h"
#include "glutil.h"
#include "vecmath.h"

#include "glm.h"

using std::min;
using std::max;
#define M_PI	3.14159265358979323846f

static bool isdrawtorus = true;
static float runtime = 0.0f;
bool specialEvent = false;
bool pause = false;
bool keyboardFlag = 0;
// ����shader���� ���ڻ��Ƴ���
GLuint drawShader;
// ���� shadowmap ���ڻ��ƹ����Ӱ
GLuint shadowmapShader; 
// �ݴ�texture
GLuint tex;
//���� ������͹�Դ ���÷��� ����ʵ����ת������ ��ʱ
void setCamera(Mtx4f &viewMatrix, Mtx4f &projectionMatrix);
void setLight(Mtx4f &viewMatrix, Mtx4f &projectionMatrix);
// ����� ����
static float camDis = 30.0f;
// �����y�� ������ת��
static float camYAng = M_PI / 4.0f;
// �������xzƽ���theta������
static float camThetaAng = M_PI / 6.0f;
// �ƹ��λ��
static Vec3f lightPos; 
//�ƹ��x����ת��
static double lightXAng = 29.9;
static double lightWaveLength = 400.0f;
//�ذ���Բ���ĵ�ǰ���� local Matrix

Mtx4f torusModelMat;

//��Ӱ��ͼ/ frame buffer/ Ĭ�Ϸֱ��� 1024
GLuint shadowmapTex; 
GLuint shadowmapFbuff; 
const int shadowmap_resolution = 4096;

GLMmodel* pmodel = NULL;


double sourcePathData[21][31][2];
double zenithAngle[21];
double waveLength[31];

int i = (lightXAng-29.9)*100;
int j = (lightWaveLength-400)/10;
float fuliangdu ;
float tongguolv ;

static void loadSourcePathData()
{
	FILE *f;
	int i=-1,j=-1,r;
	double d;
    f=fopen("Radiance Parameter for Source Path.txt","r");

    while (1) {
        r=fscanf(f,"%lf",&d);
		if (1==r){
			if (d<31 && d>29)
			{
				
				j=-1;
				i++;
				zenithAngle[i] = d;
				
				continue;
			}else
				if(d<=700 && d>=400){
					j++;
					waveLength[j]=d;
					
					continue;
				}else
				{
					if(d>1)
					sourcePathData[i][j][1] = d;
					else sourcePathData[i][j][2] = d;

					//printf("%f\n",d);
				}

		}else if (0==r) {
            fscanf(f,"%*c");
        } else break;
		
    }
fclose(f);
}

static void loadTex();

static void readPix()
{

	int x=0,y=0;
	GLfloat *Pixel = new float; 
	for(x=0;x<511;x++)
		for(y=0;y<511;y++)
		{
			glReadPixels(x,y,1,1,GL_COLOR_INDEX,GL_FLOAT,Pixel);
			printf("%f",Pixel);
		}

}
////////////////////////////
//��ʼ��
///////////////////////////
static void init()
{
	//--��ʼ�� glew �� IL
	glewInit();				
	ilInit();					
	ilutRenderer(ILUT_OPENGL);  

	glBindFragDataLocation = glBindFragDataLocationEXT;
	
	loadSourcePathData();

	//--����shader���� ����shader����
	drawShader = loadShaderProgram("shading.vert", "shading.frag");
	glBindAttribLocation(drawShader, 0, "vertex"); 	
	glBindAttribLocation(drawShader, 2, "texCoordIn");
	glBindAttribLocation(drawShader, 3, "normalIn");
	glBindFragDataLocation(drawShader, 0, "fragmentColor");
	linkShaderProgram(drawShader);

	shadowmapShader = loadShaderProgram("simple.vert", "simple.frag");
	glBindAttribLocation(shadowmapShader, 0, "vertex"); 	
	glBindFragDataLocation(shadowmapShader, 0, "fragmentColor");
	linkShaderProgram(shadowmapShader);
	
	glUseProgram( drawShader );					

	//--- ����model����
	Mtx4f ROT; 
	Mtx4f TRA; 
	Mtx4f SCA;
	SCA.scale(10.0f,10.0f,10.0f);
	ROT.rotX(M_PI/(-2.0f));
	TRA.translate(0.0f, 1.0f, 0.0f);
	torusModelMat = TRA * SCA;

	//��ʼ�� shadow texture �� frame buffer
	glGenTextures( 1, &shadowmapTex );
	glBindTexture( GL_TEXTURE_2D, shadowmapTex );		
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, shadowmap_resolution, shadowmap_resolution, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	Vec4f ones(1.0, 1.0, 1.0, 1.0); 
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, ones.vec); 

	glGenFramebuffers(1, &shadowmapFbuff);
	glBindFramebuffer(GL_FRAMEBUFFER, shadowmapFbuff);	
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadowmapTex, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE); 
	
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

}


//////////////////////
//���Ƴ���
//////////////////////
static void drawScene(const Mtx4f &viewMatrix, const Mtx4f &projectionMatrix,
					  const Mtx4f &lightViewMatrix, const Mtx4f &lightProjectionMatrix)
{
	// ���� ���Ƴ�����shader
	glUseProgram( drawShader );				

	// ����Դ��modelviewת�����󴫸�shader
	Vec3f lightpos_mv = viewMatrix.multPnt(lightPos); 
	glUniform3fv(glGetUniformLocation(drawShader, "lightPosition"), 1, lightpos_mv.vec);

	
		//����������ͼ ���ڵذ���
	int texLoc = glGetUniformLocation( drawShader, "tex0" );
	glUniform1i( texLoc, 0 );
	glActiveTexture(GL_TEXTURE0);
	
	glBindTexture(GL_TEXTURE_2D, tex);
	tex = ilutGLLoadImage("color.png");

	
		// ��shadowmap����Texture1����Ŀ����
	glUniform1i( glGetUniformLocation( drawShader, "tex1" ), 1 );
	glActiveTexture(GL_TEXTURE1); 
	glBindTexture(GL_TEXTURE_2D, shadowmapTex);

	//MV����
	Mtx4f modelViewMatrix = viewMatrix * torusModelMat;	
	Mtx4f lightMatrix = viewMatrix; 
	lightMatrix.invert(); 
	lightMatrix = lightProjectionMatrix * lightViewMatrix * lightMatrix; 
	glUniformMatrix4fv(glGetUniformLocation(drawShader, "lightMatrix"), 1, false, lightMatrix.array);
	//����MVP����
	Mtx4f modelViewProjectionMatrix = projectionMatrix * modelViewMatrix;
	Mtx4f normalMatrix = torusModelMat;
	normalMatrix.invert(); 
	normalMatrix.transpose(); 
	glUniformMatrix4fv(glGetUniformLocation(drawShader, "modelViewMatrix"),
		1, false, modelViewMatrix.array);
	glUniformMatrix4fv(glGetUniformLocation(drawShader, "modelViewProjectionMatrix"), 
		1, false, modelViewProjectionMatrix.array);
	glUniformMatrix4fv(glGetUniformLocation(drawShader, "normalMatrix"),
		1, false, normalMatrix.array);
	if (isdrawtorus)
	{
		glmDraw(pmodel, GLM_FLAT);
	}

	
	int i = (lightXAng-29.9)*100;
	int j = (lightWaveLength-400)/10;
	float fuliangdu = sourcePathData[i][j][1];
	float tongguolv = sourcePathData[i][j][2];

	float normalizeFuTong = fuliangdu * tongguolv / 51.767;
	if(keyboardFlag == 1)
	{
		printf("�춥�ǣ�%f ������%f �����ȣ�%f ͨ����%f ��һ����%f\n",lightXAng,lightWaveLength,fuliangdu,tongguolv,normalizeFuTong);
		keyboardFlag = 0;
	}

	int loc1 = glGetUniformLocation(drawShader,"fu");  
	glUniform1f(loc1,normalizeFuTong);  


	//��յ�ǰshader
	glUseProgram( 0 );	

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);

	char *c;
	char *string[11]= {"0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"};

	for (int i = 0;i<11;i++)
	{
		//glColor3f(1.0, 0.0, 0.0);
		glRasterPos2f(-0.05f+i*0.1f,0.85f);
		//glutBitmapCharacter( GLUT_BITMAP_8_BY_13, 48+i);
		for (c = string[i]; *c != '\0'; c++) 
		{
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *c);
		}

	}

	glBegin(GL_QUADS);
	glVertex2f(-0.05,0.9);
	glTexCoord2f(1,1);

	glVertex2f(0.95,0.9);
	glTexCoord2f(1,1);

	glVertex2f(0.95,1);
	glTexCoord2f(0.0,0.0);

	glVertex2f(-0.05,1);
	glTexCoord2f(0.0,0.0);
	



	glEnd();
	glEnable(GL_DEPTH_TEST);

}		


//////////////////
//������Ӱ��ͼ
/////////////////
static void drawShadowMap(const Mtx4f &viewMatrix, const Mtx4f &projectionMatrix)
{
	glPolygonOffset(2.5, 10);
	glEnable(GL_POLYGON_OFFSET_FILL);

	glBindFramebuffer(GL_FRAMEBUFFER, shadowmapFbuff);
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,shadowmap_resolution, shadowmap_resolution); 

	glClearColor(1.0, 1.0, 1.0, 1.0);//��Ӱͼ��������ɫ
	glClearDepth(1.0);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	//������ɫ��Ϊsampleshader׼��������Ӱ��ͼ
	GLint current_shader; 
	glGetIntegerv(GL_CURRENT_PROGRAM, &current_shader);
	glUseProgram( shadowmapShader );
	
	//����ģ��
	Mtx4f modelViewMatrix = viewMatrix * torusModelMat;	
	Mtx4f modelViewProjectionMatrix = projectionMatrix * modelViewMatrix;
	glUniformMatrix4fv(glGetUniformLocation(shadowmapShader, "modelViewProjectionMatrix"), 
		1, false, modelViewProjectionMatrix.array);
	//drawTorus(1.0, 0.4, 10, 40);
	glmDraw(pmodel, GLM_FLAT);
	//�л���Ըshader
	glUseProgram( current_shader );	

	glPopAttrib();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glDisable(GL_POLYGON_OFFSET_FILL);
}

/////////////////
//���
/////////////////
static void drawLight(const Mtx4f &modelMatrix, const Mtx4f &viewMatrix, const Mtx4f &projectionMatrix)
{

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glLoadIdentity();
	glTranslatef(lightPos.vec[0],lightPos.vec[1],lightPos.vec[2]);
	drawTorus(1.0, 0.4, 10, 40);
	
	glPopMatrix();


	glColor3f(1.0,0,0);
	glLineWidth(5);
	glBegin(GL_LINES);
	glVertex3f(lightPos.vec[0],lightPos.vec[1],lightPos.vec[2]);
	glVertex3f(0,0,0);
	glEnd();

}


void display(void)
{
	static float angleX=0.0,angleY=0.0;
	int rx=0,ry=0,rz=0;
	float px=0,py=0,pz=0;
	float scale=1;
	float lx=1;

	glClearColor(1.0,1.0,1.0,1.0);						// clear�󱳾�����ɫ
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // ��ɫ��z-buffer��������

	glEnable(GL_DEPTH_TEST);	// �� Z-buffering 
	glEnable(GL_CULL_FACE);		// ��ʾ������
	//viewport
	int w = glutGet((GLenum)GLUT_WINDOW_WIDTH);
	int h = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
	glViewport(0, 0, w, h);	
	//��������͹�Դ ��ñ���ľ���
	Mtx4f viewMatrix;
	Mtx4f projectionMatrix;
	setCamera(viewMatrix, projectionMatrix);	
	Mtx4f lightViewMatrix; 
	Mtx4f lightProjectionMatrix; 
	setLight(lightViewMatrix, lightProjectionMatrix);

	if (!pmodel) {
        pmodel = glmReadOBJ("obj/tree3.obj");
        if (!pmodel) exit(0);
        glmUnitize(pmodel);
        glmFacetNormals(pmodel);
        glmVertexNormals(pmodel, 90.0);

    }
	/******�����������Ӱ��ͼ***********/

	drawShadowMap(lightViewMatrix, lightProjectionMatrix); 

	Mtx4f lightModelMatrix = lightViewMatrix; 
	lightModelMatrix.invert(); 
	/*******���************************/
	//drawLight(lightModelMatrix, viewMatrix, projectionMatrix);
	
	drawScene(viewMatrix, projectionMatrix, lightViewMatrix, lightProjectionMatrix);

	glutSwapBuffers();  // ����˫����
}

/////////////////////
//���������
////////////////////
void setCamera(Mtx4f &viewMatrix, Mtx4f &projectionMatrix)
{
	Vec3f viewat = Vec3f(0,0,0); 
	Vec3f viewup = Vec3f(0,1,0); 
	Vec3f yaxis  = Vec3f(0,1,0); 
	Vec3f viewpos = Vec3f(0,0,1); 
	Mtx3f rotmat;

	// ������ת
	viewpos = Vec3f(0,0,camDis); 
	rotmat.rotAxis(yaxis, camYAng); 
	viewpos = rotmat * viewpos; 

	// ������ת
	Vec3f tmpvec = viewpos.cross(yaxis); 
	tmpvec.normalize(); 
	rotmat.rotAxis(tmpvec, camThetaAng); 
	viewpos = rotmat * viewpos; 

	
	float w = (float)glutGet((GLenum)GLUT_WINDOW_WIDTH);
	float h = (float)glutGet((GLenum)GLUT_WINDOW_HEIGHT);

	//����view �� projection����
	projectionMatrix = perspectiveMatrix(45.0f, w / h, 0.01f, 300.0f); 
	viewMatrix = lookAt(viewpos, viewat, viewup);
}


//////////////////////////////
//���ù�Դ
///////////////////////////////
void setLight(Mtx4f &viewMatrix, Mtx4f &projectionMatrix)
{
	Vec3f viewpos = Vec3f(10.0f, 20.0f, 17.6f);
	Vec3f viewat = Vec3f(0,0,0); 
	Vec3f viewup = Vec3f(0,1,0); 
	Vec3f xaxis  = Vec3f(1,0,0); 
	Vec3f yaxis  = Vec3f(0,1,0); 
	Mtx3f rotmat;

	rotmat.rotAxis(xaxis, lightXAng/180*M_PI); 
	lightPos = rotmat * viewpos; 

	projectionMatrix = perspectiveMatrix(60.0f, 1.0, 1.0f, 200.0f); 

	viewMatrix = lookAt(lightPos, viewat, viewup);
}

/////////////////////////////
//��������
////////////////////////////
void handleKeys(unsigned char key, int x, int y)
{
	switch(key)
	{

	case 27:    /* ESC �˳�*/
		exit(0); 	
	case 32:
		isdrawtorus = !isdrawtorus ;
		break;
	}
}

void handleSpecialKeys(int key, int x, int y)
{
	switch(key)
	{
	case GLUT_KEY_LEFT:
		{
		if(lightWaveLength < 401.0 )
				lightWaveLength  = 700.0;
			else
				if( lightWaveLength > 701.0)
					lightWaveLength  = 400.0;
				else
					lightWaveLength -= 10.0f ;
		//printf("%f \n",lightWaveLength); 
		keyboardFlag = 1;
		}
		break;
	case GLUT_KEY_RIGHT:
		{
		if(lightWaveLength < 390.0 )
				lightWaveLength  = 700.0;
			else
				if( lightWaveLength > 690.0)
					lightWaveLength  = 400.0;
				else
					lightWaveLength += 10.0f ;
		//printf("%f \n",lightWaveLength); 
		keyboardFlag = 1;
		};
		break;
	case GLUT_KEY_UP:
		{
			if(lightXAng<=29.89999)
				lightXAng = 30.10;
			else
				if( lightXAng>=30.09000)
					lightXAng = 29.9;
				else
				lightXAng += 0.01f;
		//printf("%f \n",lightXAng);
		keyboardFlag = 1;
		}
		break;
	case GLUT_KEY_DOWN:
		{
			if(lightXAng<=29.90 )
				lightXAng = 30.10;
			else
				if( lightXAng>30.10)
					lightXAng = 29.9;
				else
				lightXAng -= 0.01f ;

		//printf("%f \n",lightXAng);
		keyboardFlag = 1;
		}
		break;
	}
		glutPostRedisplay();
}


/////////////////////////
//������ 
////////////////////////
static bool leftDown = false;
static bool middleDown = false;
static bool rightDown = false;

static int prev_x = 0;
static int prev_y = 0;

void mouse(int button, int state, int x, int y)
{
	prev_x = x;
	prev_y = y;

	bool buttonDown = state == GLUT_DOWN;

	switch(button)
	{
	case GLUT_LEFT_BUTTON:
		if(leftDown != buttonDown)
			specialEvent = !specialEvent;
		leftDown = buttonDown;
		break;
	case GLUT_MIDDLE_BUTTON:
		middleDown = buttonDown;
		break;
	case GLUT_RIGHT_BUTTON: 
		rightDown = buttonDown;
	default:
		break;
	}
}


//////////////////////
//����˶� ʵ������ �� ��ת
//////////////////////
void motion(int x, int y)
{
	int delta_x = x - prev_x;
	int delta_y = y - prev_y;

	// ����
	if(middleDown)
	{
		camDis -= float(delta_y) * 0.3f;
		// make sure cameraDistance does not become too small
		camDis = max(0.1f, camDis);
	}
	// ��ת
	if(leftDown)
	{
		camYAng	+= float(delta_x) * 0.3f * M_PI / 180.0f;
		camThetaAng += float(delta_y) * 0.3f * M_PI / 180.0f;

		camThetaAng  = min<float>(M_PI / 2.0f - 0.05f, camThetaAng);
		camThetaAng  = max<float>(-M_PI / 2.0f + 0.05f, camThetaAng);

		while (camYAng >  M_PI * 2.0f) 
		{
			camYAng -= M_PI * 2.0f;
		}
		while (camYAng < -M_PI * 2.0f) 
		{
			camYAng += M_PI * 2.0f;
		}
	}

	prev_x = x;
	prev_y = y;
}


/////////////////////////////////
//���д���
/////////////////////////////////
void idle( void )
{
	static float startTime = float(glutGet(GLUT_ELAPSED_TIME)) / 1000.0f;
	
	if (!pause)
	{
		runtime = float(glutGet(GLUT_ELAPSED_TIME)) / 1000.0f - startTime;
	}

	//TODO:�Զ��˶��߼�
//	lightYAng = 0.05 * runtime; 

	glutPostRedisplay(); 
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	/*˫���� RGB z-buffer��ȼ��*/
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(512,512);
	glutCreateWindow("Shadow Map Lab");

	//����������ƻص�����
	glutKeyboardFunc(handleKeys);
	glutSpecialFunc(handleSpecialKeys);
	glutDisplayFunc(display);	// ��ʾ
	glutMouseFunc(mouse);		// �����
	glutMotionFunc(motion);		// ����˶�
	glutIdleFunc( idle );


	init();

	glutMainLoop();  

	

	return 0;          
}
