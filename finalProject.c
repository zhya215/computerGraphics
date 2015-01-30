#include "CSCIx229.h"
double dim=10;   //  Size of world
int th=0;         //  Azimuth of view angle
int ph=0;         //  Elevation of view angle
double asp=1;     //  Aspect ratio
int axes=1;       //  Display axes
double fov=55;       //  Field of view (for perspective)
int mode=1;
double ox = 0;
double oy = 0;
double oz = 0;
double Ex = 0;   //  Eye
double Ey = 1;   //  Eye
double Ez = 5;   //  Eye
float xrot=0; //rotate left or right
float xlrot=0; //rotate of light direction
//float yrot=0; //rotate up
double radius=0; //radius of the eye's view.

// Light values
int emission=0;
double zh=0;//rotation of another light
int Th=0,Ph=30;   //  another Light angles
//position of the light
double lx=0; //default 0
double ly=1; //default 1
double lz=3; //default 3
//spot light direction
double ldx=0;//default 0
double ldy=0;//default 0
double ldz=-1;//defualt -1
int inc       =  10;  // Ball increment
float shinyvec[]={1};
float sco=45; // spot cutoff angle

int sky[6]; //skybox texture
int bricks; //bricks texture
int floors; //floor texture
int metal; //metal texture
int candle; // candle texture
int stop=0; //stop the animation

double speed=0.1; //speed of the light and perspective

double routine[9][2]={{0, 1}, {8,1}, {8,-3}, {2,-3}, {2,-7}, {-4,-7}, {-4,-11}, {0,-11}, {0,-14}};
double positionDroplight[6][2]={{0,5}, {4,1}, {8,-1}, {4,-3}, {-2,-7}, {-2,-11}};


int num_points_eye=0;//number of the points the eyes have reached.
int num_points_lig=0;//number of the points the light has reached.

double rot1=0; //for the rotation of droplight by axis z
double rot2=90;//for the rotation of droplight by axis x

int switchScene=0;//switch to the other secene
int switchObject=0;

GLuint candlestick0; //displaylist for a candlestick
GLuint candlestick1; //displaylist for a candlestick

//for the direction of eyes and light
#define TURN_LEFT -90 
#define TURN_RIGHT 90 


void rotateAroundY(double angle){
    double mat[16];
    mat[0] = Cos(angle);   mat[4] = 0;   mat[ 8] = Sin(angle);   mat[12] = 0;
    mat[1] = 0;   mat[5] = 1;    mat[ 9] = 0;   mat[13] = 0;
    mat[2] = -Sin(angle);            mat[6] = 0;             mat[10] = Cos(angle);   mat[14] = 0;
    mat[3] = 0;            mat[7] = 0;             mat[11] = 0;   mat[15] = 1;
    
    glMultMatrixd(mat);
}

void rotateAroundX(double angle){
    double mat[16];
    mat[0] = 1;   mat[4] = 0;   mat[ 8] = 0;   mat[12] = 0;
    mat[1] = 0;   mat[5] = Cos(angle);   mat[ 9] = -Sin(angle);   mat[13] = 0;
    mat[2] = 0;   mat[6] = Sin(angle);   mat[10] = Cos(angle);   mat[14] = 0;
    mat[3] =  0;   mat[7] =  0;   mat[11] =  0;   mat[15] = 1;
    
    glMultMatrixd(mat);
    
}
void rotateAroundZ(double angle){
    double mat[16];
    mat[0] = Cos(angle);   mat[4] = -Sin(angle);   mat[ 8] = 0;   mat[12] = 0;
    mat[1] = Sin(angle);   mat[5] = Cos(angle);    mat[ 9] = 0;   mat[13] = 0;
    mat[2] = 0;            mat[6] = 0;             mat[10] = 1;   mat[14] = 0;
    mat[3] = 0;            mat[7] = 0;             mat[11] = 0;   mat[15] = 1;
    glMultMatrixd(mat);
}
/*
 *  Draw a ball
 *     at (x,y,z)
 *     radius r
 */
static void ball(double x,double y,double z,double r)
{
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  Yellow ball
   glColor3f(1,1,0);
   glutSolidSphere(1.0,16,16);
   //  Undo transofrmations
   glPopMatrix();
}

/*
 *  Draw vertex in polar coordinates with normal
 */
static void Vertex(double th,double ph)
{
    double x = Sin(th)*Cos(ph);
    double y = Cos(th)*Cos(ph);
    double z =         Sin(ph);
    //  For a sphere at the origin, the position
    //  and normal vectors are the same
    glNormal3d(x,y,z);
    glVertex3d(x,y,z);
}

/*Draw a ball at (x,y,z)
 *with radius "r"
 *can be rotated by axis x
 *can be limited by the ph angle
 *can change into differen textures
 *can choose wheter to be emission
 */
static void Ball(double x,double y,double z,double r, double angle, double ph_angle, int texture, int emission)
{
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
    int th,ph;
    float yellow[] = {1.0,1.0,0.0,1.0};
    double emission_candle=0;
    if (emission == 1){
        emission_candle=100;
    }
    float Emission[]  = {0.01*emission_candle,0.0,0.0,1.0};
    
    glPushMatrix();
    glTranslated(x, y, z);
    rotateAroundX(angle);
    //  Offset, scale and rotate
    glScaled(r,r,r);
    //  White ball
    glColor3f(1,1,1);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
    
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D,texture);
    double t_width_step=1.0/(60.0/inc);
    double t_height_step=1.0/(360.0/(2*inc));
    double t_width=0;
    //  Bands of latitude
    for (ph=-90;ph<-90+ph_angle;ph+=inc)
    {
        double t_height=0;
        glBegin(GL_QUAD_STRIP);
        for (th=0;th<=360;th+=2*inc)
        {
            glTexCoord2f(t_width,t_height);
            Vertex(th,ph);
            glTexCoord2f(t_width,t_height+t_height_step);
            Vertex(th,ph+inc);
            t_height+=t_height_step;
        }
        glEnd();
        t_width+=t_width_step;
    }
    glDisable(GL_TEXTURE_2D);
    //  Undo transofrmations
    glPopMatrix();
}


static void Sky(double D)
{
   glColor3f(1,1,1);
   glEnable(GL_TEXTURE_2D);

   //  Sides
   glBindTexture(GL_TEXTURE_2D, sky[5]);
   glBegin(GL_QUADS);
   glTexCoord2f(0.00,0); glVertex3f(-D,-D,-D);
   glTexCoord2f(1,0); glVertex3f(+D,-D,-D);
   glTexCoord2f(1,1); glVertex3f(+D,+D,-D);
   glTexCoord2f(0,1); glVertex3f(-D,+D,-D);
   glDisable(GL_TEXTURE_2D);
   glEnd();

   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, sky[4]);
   glBegin(GL_QUADS);
   glTexCoord2f(0,0); glVertex3f(+D,-D,-D);
   glTexCoord2f(1,0); glVertex3f(+D,-D,+D);
   glTexCoord2f(1,1); glVertex3f(+D,+D,+D);
   glTexCoord2f(0,1); glVertex3f(+D,+D,-D);
   glDisable(GL_TEXTURE_2D);
   glEnd();

   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, sky[3]);
   glBegin(GL_QUADS);
   glTexCoord2f(0,0); glVertex3f(+D,-D,+D);
   glTexCoord2f(1,0); glVertex3f(-D,-D,+D);
   glTexCoord2f(1,1); glVertex3f(-D,+D,+D);
   glTexCoord2f(0,1); glVertex3f(+D,+D,+D);
   glDisable(GL_TEXTURE_2D);
   glEnd();

   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, sky[2]);
   glBegin(GL_QUADS);   
   glTexCoord2f(0,0); glVertex3f(-D,-D,+D);
   glTexCoord2f(1.00,0); glVertex3f(-D,-D,-D);
   glTexCoord2f(1.00,1); glVertex3f(-D,+D,-D);
   glTexCoord2f(0,1); glVertex3f(-D,+D,+D);
   glDisable(GL_TEXTURE_2D);
   glEnd();

   //  Top and bottom
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D,sky[0]);
   glBegin(GL_QUADS);
   glTexCoord2f(0,0); glVertex3f(+D,+D,-D);
   glTexCoord2f(1,0); glVertex3f(+D,+D,+D);
   glTexCoord2f(1,1); glVertex3f(-D,+D,+D);
   glTexCoord2f(0,1); glVertex3f(-D,+D,-D);
   glDisable(GL_TEXTURE_2D);
   glEnd();

   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D, sky[1]);
   glBegin(GL_QUADS);
   glTexCoord2f(0,0); glVertex3f(-D,-D,+D);
   glTexCoord2f(1,0); glVertex3f(+D,-D,+D);
   glTexCoord2f(1,1); glVertex3f(+D,-D,-D);
   glTexCoord2f(0,1); glVertex3f(-D,-D,-D);
   glEnd();

   glDisable(GL_TEXTURE_2D);
}



/*draw a circular ring at postion(x, y, z)
 *with radius R
 *and small radius r
 */
static void circularRing(double x, double y, double z, double R, double r){
    int horizontalAngle=0;
    int verticalAngle=0;
    int d=5;
    float yellow[] = {1.0,1.0,0.0,1.0};
    float Emission[]  = {0.0,0.0,0.01*emission,1.0};
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    //  Save transformation
    glPushMatrix();
    //  Offset, scale and rotate
    glTranslated(x,y,z);
    //glScaled(R+r,r,R+r);
    
    glColor3f(1,1,1);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
    
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D,metal);
    
    double t_width_step=1.0/(360/d);
    double t_height_step=1.0/(360/d);
    double t_w=0;
    
    for(verticalAngle=0; verticalAngle <= 360; verticalAngle+=d){
        glBegin(GL_QUADS);
        double t_h=0;
        for(horizontalAngle=0; horizontalAngle <= 360; horizontalAngle+=d){
            //glBegin(GL_QUAD_STRIP);
            double centreX0=Cos(horizontalAngle)*R;
            double centreZ0=Sin(horizontalAngle)*R;
            
            double centreX1=Cos(horizontalAngle+d)*R;
            double centreZ1=Sin(horizontalAngle+d)*R;
            
            double x00=centreX0-r*Cos(verticalAngle)*Cos(horizontalAngle);
            double y00=r*Sin(verticalAngle);
            double z00=centreZ0-r*Cos(verticalAngle)*Sin(horizontalAngle);
            
            double x01=centreX0-r*Cos(verticalAngle+d)*Cos(horizontalAngle);
            double y01=r*Sin(verticalAngle+d);
            double z01=centreZ0-r*Cos(verticalAngle+d)*Sin(horizontalAngle);
            
            double x10=centreX1-r*Cos(verticalAngle)*Cos(horizontalAngle+d);
            double y10=r*Sin(verticalAngle);
            double z10=centreZ1-r*Cos(verticalAngle)*Sin(horizontalAngle+d);
            
            double x11=centreX1-r*Cos(verticalAngle+d)*Cos(horizontalAngle+d);
            double y11=r*Sin(verticalAngle+d);
            double z11=centreZ1-r*Cos(verticalAngle+d)*Sin(horizontalAngle+d);
            
            glNormal3d(x00-centreX0, y00, z00-centreZ0);
            glTexCoord2f(t_w, t_h);
            glVertex3d(x00, y00, z00);
            
            glNormal3d(x01-centreX0, y01, z01-centreZ0);
            glTexCoord2f(t_w+t_width_step, t_h);
            glVertex3d(x01, y01, z01);
            
            glNormal3d(x11-centreX1, y11, z11-centreZ1);
            glTexCoord2f(t_w+t_width_step, t_h+t_height_step);
            glVertex3d(x11, y11, z11);
            
            glNormal3d(x10-centreX1, y10, z10-centreZ1);
            glTexCoord2f(t_w, t_h+t_height_step);
            glVertex3d(x10, y10, z10);
            
            t_h+=t_height_step;
        }
        glEnd();
    }
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

/*r1 begin radius r2 end radius
 *rorate by x axis
 *horizontal angle angleH
 */
static void halfRing(double x, double y, double z, double R, double r1, double r2,double dx, double dy, double dz, double angle, int angleH){
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    int horizontalAngle=0;
    int verticalAngle=0;
    int d=5;
    float yellow[] = {1.0,1.0,0.0,1.0};
    float Emission[]  = {0.0,0.0,0.01*emission,1.0};
    double step=(r1-r2)/(angleH/d);
    glPushMatrix();
    glTranslated(x, y, z);
    rotateAroundX(angle);
    
    //  Offset, scale and rotate
    glScaled(dx,dy,dz);
    
    glColor3f(1,1,1);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
    
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D,metal);
    
    double t_width_step=1.0/(angleH/d);
    double t_height_step=1.0/(360.0/d);
    
    double t_w=0;
    
    for(verticalAngle=0; verticalAngle <= 360; verticalAngle+=d){
        glBegin(GL_QUADS);
        double r=r2;
        double t_h=0;
        t_w+=t_width_step;
        for(horizontalAngle=0; horizontalAngle <= angleH; horizontalAngle+=d){
            //glBegin(GL_QUAD_STRIP);
            //r+=step;
            double centreX0=Cos(horizontalAngle)*R;
            double centreZ0=Sin(horizontalAngle)*R;
            
            double centreX1=Cos(horizontalAngle+d)*R;
            double centreZ1=Sin(horizontalAngle+d)*R;
            
            double x00=centreX0-r*Cos(verticalAngle)*Cos(horizontalAngle);
            double y00=r*Sin(verticalAngle);
            double z00=centreZ0-r*Cos(verticalAngle)*Sin(horizontalAngle);
            
            double x01=centreX0-r*Cos(verticalAngle+d)*Cos(horizontalAngle);
            double y01=r*Sin(verticalAngle+d);
            double z01=centreZ0-r*Cos(verticalAngle+d)*Sin(horizontalAngle);
            
            double x10=centreX1-(r+step)*Cos(verticalAngle)*Cos(horizontalAngle+d);
            double y10=(r+step)*Sin(verticalAngle);
            double z10=centreZ1-(r+step)*Cos(verticalAngle)*Sin(horizontalAngle+d);
            
            double x11=centreX1-(r+step)*Cos(verticalAngle+d)*Cos(horizontalAngle+d);
            double y11=(r+step)*Sin(verticalAngle+d);
            double z11=centreZ1-(r+step)*Cos(verticalAngle+d)*Sin(horizontalAngle+d);
            
            r+=step;
            
            glNormal3d(x00-centreX0, y00, z00-centreZ0);
            glTexCoord2f(t_w, t_h);
            glVertex3d(x00, y00, z00);
            
            glNormal3d(x01-centreX0, y01, z01-centreZ0);
            glTexCoord2f(t_w+t_width_step, t_h);
            glVertex3d(x01, y01, z01);
            
            glNormal3d(x11-centreX1, y11, z11-centreZ1);
            glTexCoord2f(t_w+t_width_step, t_h+t_height_step);
            glVertex3d(x11, y11, z11);
            
            glNormal3d(x10-centreX1, y10, z10-centreZ1);
            glTexCoord2f(t_w, t_h+t_height_step);
            glVertex3d(x10, y10, z10);
            t_h+=t_height_step;
        }
        glEnd();
    }
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}
/*
 *Rotate halfring by z axis
 */
static void rotateRingZ(double x, double y, double z, double angle){
    
    glPushMatrix();
    glTranslated(x, y, z);
    rotateAroundZ(angle);
    
    halfRing(0,0,0,1,0.1,0, 1, 1, 1, 90, 180);
    
    glPopMatrix();
}
/*
 *Rotate halfring by y axis
 */
static void rotateRingY(double x, double y, double z, double angle){
    
    glPushMatrix();
    glTranslated(x, y, z);
    rotateAroundY(angle);
    
    
    halfRing(0,0,0,1.5,0.3,0.1, 1, 1, 1, 90, 180);
    
    glPopMatrix();
}

/*
 *Rotate circularring by z axis
 */
static void rotateCircZ(double x, double y, double z, double R, double r, double angle){
    glPushMatrix();
    glTranslated(x,y,z);
    rotateAroundZ(angle);
    
    circularRing(0,0,0,R,r);
    glPopMatrix();
    
}

/*
 *Rotate circularring by x axis
 */
static void rotateCircX(double x, double y, double z, double R, double r, double angle){
    glPushMatrix();
    glTranslated(x,y,z);
    rotateAroundX(angle);
    
    circularRing(0,0,0,R,r);
    glPopMatrix();
    //glEnd();
}

/*Draw a hander for candlestick at position (x, y, z)
 *with the scale (dx, dy, dz)
 *rotate by axis y
 */
void hander(double x, double y, double z, double dx, double dy, double dz, double angle)
{
    
    glPushMatrix();
    glTranslated(x, y, z);
    glScaled(dx, dy, dz);
    rotateAroundY(angle);
    
    halfRing(1,0,0,1,0.2,0.1, 1, 1, 1, -90, 200);
    rotateRingZ(3,0,0, 0);
    
    glPopMatrix();
}

/*
 *Draw a the centre hander at postion (x,y,z)
 *with the top radius is r and buttom is R
 *with the height h
 *can be rotate by z axis
 */
void centreHander(double x, double y, double z, double r, double R, double h, double angle){
    float yellow[] = {1.0,1.0,0.0,1.0};
    float Emission[]  = {0.0,0.0,0.01*emission,1.0};
    int step=5;
    int i=0;
    
    glPushMatrix();
    glTranslated(x, y, z);
    rotateAroundZ(angle);
    
    glColor3f(1,1,1);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
    
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D,metal);
    
    
    double t_width_step=1.0/(360/step);
    double t_width=0;
    
    glBegin(GL_QUADS);
    for(i=0; i<360; i+=step){
        double xr0=r*Cos(i);
        double yr0=0;
        double zr0=r*Sin(i);
        
        double xr1=r*Cos(i+step);
        double yr1=0;
        double zr1=r*Sin(i+step);
        
        double xR0=R*Cos(i);
        double yR0=-h;
        double zR0=R*Sin(i);
        
        double xR1=R*Cos(i+step);
        double yR1=-h;
        double zR1=R*Sin(i+step);
        
        
        glNormal3d(h*Cos(i), (R-r), h*Sin(i));
        glTexCoord2f(0, t_width);
        glVertex3d(xr0, yr0, zr0);
        
        
        glNormal3d(h*Cos(i+step), (R-r), h*Sin(i+step));
        glTexCoord2f(0, t_width+t_width_step);
        glVertex3d(xr1, yr1, zr1);
        
        
        glNormal3d(h*Cos(i+step), (R-r), h*Sin(i+step));
        glTexCoord2f( 1, t_width+t_width_step);
        glVertex3d(xR1, yR1, zR1);
        
        
        glNormal3d(h*Cos(i), (R-r), h*Sin(i));
        glTexCoord2f(1, t_width);
        glVertex3d(xR0, yR0, zR0);
        
        t_width+=t_width_step;
    }
    
    glEnd();
    glDisable(GL_TEXTURE_2D);
    
    
    glPopMatrix();
}

/*Draw a cylinder
 *
 */
void cylinder(double x, double y, double z, int texture, double r, double h){
    glPushMatrix();
    glTranslated(x,y,z);
    int step=5;
    int i=0;
    float Emission[]  = {0.0,0.0,0.0,1.0};
    glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
    
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
    glBindTexture(GL_TEXTURE_2D,texture);
    
    double t_width_step=1.0/(360/step);
    double t_width=0;
    
    glBegin(GL_QUADS);
    for(i=0; i<=360; i+=step){
        double xd0=r*Cos(i);
        double yd0=0;
        double zd0=r*Sin(i);
        
        double xd1=r*Cos(i+step);
        double yd1=0;
        double zd1=r*Sin(i+step);
        
        double xu0=r*Cos(i);
        double yu0=h;
        double zu0=r*Sin(i);
        
        double xu1=r*Cos(i+step);
        double yu1=h;
        double zu1=r*Sin(i+step);
        
        glNormal3d(xd0, yd0, zd0);
        glTexCoord2f(t_width, 0);
        glVertex3d(xd0, yd0, zd0);
        
        glNormal3d(xd1, yd1, zd1);
        glTexCoord2f(t_width+t_width_step, 0);
        glVertex3d(xd1, yd1, zd1);
        
        glNormal3d(xd1, yd1, zd1);
        glTexCoord2f(t_width+t_width_step, 1);
        glVertex3d(xu1, yu1, zu1);
        
        glNormal3d(xd0, yd0, zd0);
        glTexCoord2f(t_width, 1);
        glVertex3d(xu0, yu0, zu0);
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

/*Draw a candle at position (x,y,z)
 *with height h
 *radius r
 */
void Candle(double x, double y, double z, double h, double r){
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    float yellow[] = {1.0,1.0,0.0,1.0};
    float Emission[]  = {0.0,0.0,0.01*emission,1.0};
    
    
    glPushMatrix();
    glTranslated(x, y, z);
    
    glColor3f(1,1,1);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
    cylinder(0,0,0,candle,r,h);
    Ball(0,h+r/sqrt(3),0,r*2/sqrt(3), -90, 60, candle, 0);
    Ball(0, h+r/sqrt(3),0,r*0.4, 0, 180, candle, 1);
    
    glPopMatrix();
}

/* Draw a part of chain at position (x,y,z)
 *with high "h", width "w" and radius "r"
 *can be rotated by y axis
 */
void partChain(double x, double y, double z, double r, double h, double w, double angle){
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    float yellow[] = {1.0,1.0,0.0,1.0};
    
    glPushMatrix();
    glTranslated(x,y,z);
    rotateAroundY(angle);
    glColor3f(1,1,1);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    cylinder(-w/2,-h/2,0,metal,r,h);
    cylinder(w/2, -h/2, 0, metal, r, h);
    halfRing(0,h/2,0,w/2,r,r, 1, 1, 1, -90, 180);
    halfRing(0,-h/2,0,w/2,r,r, 1, 1, 1, 90, 180);
    glPopMatrix();
}

/*Draw a candlestick at position (x,y,z)
 *scale by (dx, dy, dz)
 *can be rotate by axis y
 */
static void candlestick(double x, double y, double z, double dx, double dy, double dz, double angle){
    
    glPushMatrix();
    glTranslated(x, y, z);
    glScaled(dx, dy, dz);
    rotateAroundY(angle);
    
    
    Ball(0,2.8,0, 1.5,-90, 60, metal, 0);
    Ball(-4.5,1.3,0, 1.5, -90, 60, metal, 0);
    Ball(4.5,1.3,0, 1.5, -90, 60, metal, 0);
    circularRing(0,-1,0,0.3,0.1);
    hander(0.6,0,0,1,1,1,0);
    hander(-0.6,0,0,1,1,1,180);
    //centreHander(0,-1,0,0.3,1,4,0);
    centreHander(0,-1,0,0.3,0.7,2.5,180);
    Candle(0,1.5,0,1,0.5);
    Candle(-4.5, 0, 0, 1, 0.5);
    Candle(4.5, 0, 0, 1, 0.5);
    rotateRingY(0,-1.2,-1.5,90);
    
    glPopMatrix();
}

/*Draw a chain at position (x,y,z)
 *
 */
static void chain(double x, double y, double z){
    glPushMatrix();
    glTranslated(x,y,z);
    partChain(0,0,0, 0.1, 0.5, 0.5, 0);
    partChain(0,0.75,0, 0.1, 0.5, 0.5, 90);
    partChain(0,1.5,0, 0.1, 0.5, 0.5, 0);
    partChain(0,2.25,0, 0.1, 0.5, 0.5, 90);
    partChain(0,3,0, 0.1, 0.5, 0.5, 0);
    glPopMatrix();
}

/*Draw a drop light at position (x,y,z)
 * can be scaled by (dx,dy,dz)
 *whether to turn on the light
 */
static void droplight(double x, double y, double z, double dx, double dy, double dz, int turnOn){
    glPushMatrix();
    glTranslated(x,y,z);
    glScaled(dx,dy,dz);
    
    if(turnOn){
        //Create light for the droplight
     float ambient[]={0 ,0 ,0 ,1.0};
        float diffuse[]={1 ,1 ,1 ,1.0};
        float specular[]={0,0,0,1.0};
    
        float position[]={0,0,0,1};
        float direction[]={0,-1,0,0};
    
        float lsco=45;
    
        //  Enable lighting
        glEnable(GL_LIGHTING);
        //  Location of viewer for specular calculations
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,1);
        glEnable(GL_NORMALIZE);
        glEnable(GL_LIGHT2);
    
        //  Set ambient, diffuse, specular components and position of light 0
        glLightfv(GL_LIGHT2,GL_AMBIENT ,ambient);
        glLightfv(GL_LIGHT2,GL_DIFFUSE ,diffuse);
        glLightfv(GL_LIGHT2,GL_SPECULAR,specular);
        glLightfv(GL_LIGHT2,GL_POSITION,position);
        //  Set spotlight parameters
        glLightfv(GL_LIGHT2,GL_SPOT_DIRECTION,direction);
        glLightf(GL_LIGHT2,GL_SPOT_CUTOFF,lsco);
        glLightf(GL_LIGHT2,GL_SPOT_EXPONENT,0);
    
        //  Set attenuation
        glLightf(GL_LIGHT2,GL_CONSTANT_ATTENUATION ,1);
        glLightf(GL_LIGHT2,GL_LINEAR_ATTENUATION   ,0);
        glLightf(GL_LIGHT2,GL_QUADRATIC_ATTENUATION,0);
    }
    
    rotateCircZ(0,0,0,4.5,0.2,0);
    rotateCircZ(0,0,0,4.5,0.2,90);
    
    rotateCircZ(0,-3,0,3.2,0.2,0);
    rotateCircZ(0, -2, 0, 4, 0.2, 0);
    rotateCircZ(0,-1,0,4.2,0.2,0);
    
    
    rotateCircX(0,0,0,4.5,0.2,90);
    rotateCircZ(0,0,0,2,0.1,rot1);
    rotateCircX(0,0,0,2,0.1,rot2);
    
    Ball(0,0,0, 0.5,0, 180, candle, 1);
    chain(0,5,0);
    
    
    
    glPopMatrix();
}


/*Draw a wall at position (x, y, z)
 *can be rotated by axis z
 */
static void wall(double x, double y, double z, double angle){
  double length=2.0;
  double width=2.0;
  int num=100; //divide the side of quad into 100 parts.
   int i=0;
   int j=0;

   glPushMatrix();
   glTranslated(x, y, z);
    rotateAroundZ(angle);

   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,bricks);

   //Draw a roof

   glNormal3f(0,1,0); 
   glBegin(GL_QUADS);
   double l_tex=length/num;
   double l=-2.0/num;
   double w_tex=2.0/num;
   double w=width/num;
   for(i=0; i<num; i++){
    for(j=0; j<num; j++){
      glTexCoord2f(l_tex*j, w_tex*i);
      glVertex3d(l*j, 0, w*i);
      glTexCoord2f(l_tex*j, w_tex*(i+1));
      glVertex3d(l*j, 0, w*(i+1));
      glTexCoord2f(l_tex*(j+1), w_tex*(i+1));
      glVertex3d(l*(j+1), 0, w*(i+1));
      glTexCoord2f(l_tex*(j+1), w_tex*i);
      glVertex3d(l*(j+1), 0, w*i);
    }
   }


   glEnd();

   
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();

}

/*Draw a wall at position (x, y, z)
 *can be rotated by axis x
 */
static void wallX(double x, double y, double z, double angle){
  double length=2.0;
  double width=2.0;
  int num=100; //divide the side of quad into 100 parts.
   int i=0;
   int j=0;

   glPushMatrix();
   glTranslated(x, y, z);
    rotateAroundX(angle);

   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,bricks);

   //Draw a roof

   glNormal3f(0,1,0); 
   glBegin(GL_QUADS);
   double l_tex=length/num;
   double l=-2.0/num;
   double w_tex=2.0/num;
   double w=width/num;
   for(i=0; i<num; i++){
    for(j=0; j<num; j++){
      glTexCoord2f(l_tex*j, w_tex*i);
      glVertex3d(w*i, 0, l*j);
      glTexCoord2f(l_tex*j, w_tex*(i+1));
      glVertex3d(w*(i+1), 0, l*j);
      glTexCoord2f(l_tex*(j+1), w_tex*(i+1));
      glVertex3d(w*(i+1), 0, l*(j+1));
      glTexCoord2f(l_tex*(j+1), w_tex*i);
      glVertex3d(w*i, 0, l*(j+1));
    }
   }


   glEnd();

   
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();

}

/*
 *Draw a floor at position(x,y,z)
 */
static void Floor(double x, double y, double z){
   double length=2.0;
   double width=2.0;
   int num=100; //divide the side of quad into 100 parts.
   int i=0;
   int j=0;
   glPushMatrix();
   glTranslated(x, y, z);

   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,floors);


   glNormal3f(0,1,0); 
   glBegin(GL_QUADS);
   double l_tex=length/num;
   double l=-2.0/num;
   double w_tex=2.0/num;
   double w=width/num;
   for(i=0; i<num; i++){
    for(j=0; j<num; j++){
      glTexCoord2f(l_tex*j, w_tex*i);
      glVertex3d(l*j, 0, w*i);
      glTexCoord2f(l_tex*j, w_tex*(i+1));
      glVertex3d(l*j, 0, w*(i+1));
      glTexCoord2f(l_tex*(j+1), w_tex*(i+1));
      glVertex3d(l*(j+1), 0, w*(i+1));
      glTexCoord2f(l_tex*(j+1), w_tex*i);
      glVertex3d(l*(j+1), 0, w*i);
    }
   }
   glEnd();

   glPopMatrix();
}
/* Draw top and sides of the maze at postion (x,y,z)
 * can be rotated by y axis
 * whether to show the droplight
 * whether turn on the light
 */
static void maze(double x, double y, double z, double angle, int dropl, int turnOn){

   glPushMatrix();
   glTranslated(x, y, z);
    rotateAroundY(angle);
    if(dropl){
        if(turnOn){
            droplight(0,1.6,1,0.05,0.05,0.05,1);
        }
        else{
            droplight(0,1.6,1,0.05,0.05,0.05,0);
        }
    }

    glCallList(candlestick0);
    glCallList(candlestick1);
   wall(-1,0,0,-90);
   wall(-1,2,0, -180);
   wall(1, 2, 0, -270);

   double length=2.0;
   double width=2.0;
   int num=100; //divide the side of quad into 100 parts.
   int i=0;
   int j=0;

   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,floors);


   glNormal3f(0,1,0); 
   glBegin(GL_QUADS);
   double l_tex=length/num;
   double l=-2.0/num;
   double w_tex=2.0/num;
   double w=width/num;
   for(i=0; i<num; i++){
    for(j=0; j<num; j++){
      glTexCoord2f(l_tex*j, w_tex*i);
      glVertex3d(l*j+1, 0, w*i);
      glTexCoord2f(l_tex*j, w_tex*(i+1));
      glVertex3d(l*j+1, 0, w*(i+1));
      glTexCoord2f(l_tex*(j+1), w_tex*(i+1));
      glVertex3d(l*(j+1)+1, 0, w*(i+1));
      glTexCoord2f(l_tex*(j+1), w_tex*i);
      glVertex3d(l*(j+1)+1, 0, w*i);
    }
   }
   glEnd();
   // glDisable(GL_LIGHTING);

   glPopMatrix();


}


/*
 *Init gl fog before using it
 */
static void setFog(){
    GLfloat fogColor[4] = {0.5, 0.5, 0.5, 1.0};
    GLfloat density = 0.3;
    //glEnable (GL_DEPTH_TEST);
    
    glEnable(GL_FOG);
    
    glFogi(GL_FOG_MODE, GL_EXP2);
    
    glFogfv(GL_FOG_COLOR, fogColor);
    glFogf(GL_FOG_DENSITY, density);
    glHint(GL_FOG_HINT, GL_NICEST);
}
/*
 *determine whether the eyes is in the space
 */
int ifIn(double x, double z){
    return Ex >= x-2 && Ex <= x+2 && Ez >= z-2 && Ez <= z+2;
}
/*
 *Draw the whole Maze.
 */
static void wholeMaze(){
   maze(0,0,2,0,0,0);

   maze(0,0,4,0,1,ifIn(positionDroplight[0][0], positionDroplight[0][1]));
   
   wall(-3,0,0,-90);
   //first crossing
   Floor(1,0,0);
   wallX(-1,0,0,90);
   wallX(-1,2,0,180);
    maze(1,0,1,90,0,0);

   maze(3,0,1,90,1,ifIn(positionDroplight[1][0], positionDroplight[1][1]));//right
   maze(-1,0,1,-90,0,0);//left

   maze(5,0,1,90,0,0);
   //second crossing
   Floor(9,0,0);
   wallX(7,2,0,180);
   wallX(7,2,2,270);
   wall(9,2,0,-270);
   maze(8,0,-2,0,1,ifIn(positionDroplight[2][0], positionDroplight[2][1]));//left
   //third crossing
   Floor(9,0,-4);
   wallX(7,0,-4,90);
   wallX(7,2,-4,180);
   maze(7,0,-3, -90,0,0);//left
   maze(11,0,-3, -90,0,0);//right
   wall(11,2,-4,-270);

   maze(5,0,-3,-90,1,ifIn(positionDroplight[3][0], positionDroplight[3][1]));

   //fourth crossing
   Floor(3, 0, -4);
   wallX(1, 2, -4, 180);
   wallX(1, 2, -2, 270);
   wall(1,0,-4, -90);
   maze(2,0,-6, 0,0,0);//right
   //fifth crossing
   Floor(3,0,-8);
   wallX(1, 2, -8, 180);
   wall(3, 2, -8, -270);
   maze(1,0,-7,-90,0,0);//left
   wallX(1,0,-10,90);
   maze(2,0,-10,0,0,0);//straight

   maze(-1,0,-7,-90,1,ifIn(positionDroplight[4][0], positionDroplight[4][1]));

   //sixth crossing
   Floor(-3,0,-8);
   wallX(-5, 2, -8, 180);
   wallX(-5,2,-6,270);
   wall(-5,0,-8,-90);
   maze(-4,0,-10,0,0,0);//right
   //seventh crossing
   Floor(-3,0,-12);
   wallX(-5,2,-12,180);
   wall(-5,0,-12,-90);
   wallX(-5,0,-12,90);
   maze(-3, 0, -11, 90,1, ifIn(positionDroplight[5][0],positionDroplight[5][1]));//right
   //eighth crossing
   Floor(1,0,-12);
   wallX(-1,2,-12,180);
   wallX(-1,2,-10,270);
   wall(1,2,-12,-270);

   maze(0, 0, -14, 0,0,0);
}
/*Set the light number, position and direction
 *
 */
static void setLight(int light_num, double x, double y, double z, double dx, double dy, double dz, double lightsco){
    glShadeModel(GL_SMOOTH);
    // set light
    float Ambient[]={0 ,0 ,0 ,1.0};
    float Diffuse[]={1 ,1 ,1 ,1.0};
    float Specular[]={0.05,0.05,0.05,1.0};
    
    //set color and direction of spotlight
    float yellow[] = {1.0,1.0,0.0,1.0};
    float Direction[] = {dx,dy,dz,0};
    
    //set the position of light
    float Position[]={x, y, z, 1};
    ball(Position[0], Position[1], Position[2], 0.05);
    
    //  OpenGL should normalize normal vectors
    glEnable(GL_NORMALIZE);
    //  Enable lighting
    glEnable(GL_LIGHTING);
    //  Location of viewer for specular calculations
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,1);
    //  Two sided mode
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    //  glColor sets ambient and diffuse color materials
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    //  Set specular colors
    glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
    glMaterialfv(GL_FRONT,GL_SHININESS,shinyvec);
    //  Enable light 0
    glEnable(light_num);
    //  Set ambient, diffuse, specular components and position of light 0
    glLightfv(light_num,GL_AMBIENT ,Ambient);
    glLightfv(light_num,GL_DIFFUSE ,Diffuse);
    glLightfv(light_num,GL_SPECULAR,Specular);
    glLightfv(light_num,GL_POSITION,Position);
    //  Set spotlight parameters
    glLightfv(light_num,GL_SPOT_DIRECTION,Direction);
    glLightf(light_num,GL_SPOT_CUTOFF,lightsco);
    glLightf(light_num,GL_SPOT_EXPONENT,0);
    //  Set attenuation
    glLightf(light_num,GL_CONSTANT_ATTENUATION ,1);
    glLightf(light_num,GL_LINEAR_ATTENUATION   ,0);
    glLightf(light_num,GL_QUADRATIC_ATTENUATION,0);
}

static void anotherScene(){
    switch (switchObject) {
        case 0:
            hander(0.6,0,0,1,1,1,0);
            hander(-0.6,0,0,1,1,1,180);
            break;
        case 1:
            
            hander(0.6,0,0,1,1,1,0);
            hander(-0.6,0,0,1,1,1,180);
            centreHander(0,-1,0,0.3,0.7,2.5,180);
            
            
            break;
        case 2:
            hander(0.6,0,0,1,1,1,0);
            hander(-0.6,0,0,1,1,1,180);
            centreHander(0,-1,0,0.3,0.7,2.5,180);
            circularRing(0,-1,0,0.3,0.1);
            
            break;
        case 3:
            hander(0.6,0,0,1,1,1,0);
            hander(-0.6,0,0,1,1,1,180);
            centreHander(0,-1,0,0.3,0.7,2.5,180);
            circularRing(0,-1,0,0.3,0.1);
            rotateRingY(0,-1.2,-1.5,90);
            
            
            break;
        case 4:
            hander(0.6,0,0,1,1,1,0);
            hander(-0.6,0,0,1,1,1,180);
            centreHander(0,-1,0,0.3,0.7,2.5,180);
            circularRing(0,-1,0,0.3,0.1);
            rotateRingY(0,-1.2,-1.5,90);
            Ball(0,2.8,0, 1.5,-90, 60, metal, 0);
            Ball(-4.5,1.3,0, 1.5, -90, 60, metal, 0);
            Ball(4.5,1.3,0, 1.5, -90, 60, metal, 0);
            break;
        case 5:
            hander(0.6,0,0,1,1,1,0);
            hander(-0.6,0,0,1,1,1,180);
            centreHander(0,-1,0,0.3,0.7,2.5,180);
            circularRing(0,-1,0,0.3,0.1);
            rotateRingY(0,-1.2,-1.5,90);
            Ball(0,2.8,0, 1.5,-90, 60, metal, 0);
            Ball(-4.5,1.3,0, 1.5, -90, 60, metal, 0);
            Ball(4.5,1.3,0, 1.5, -90, 60, metal, 0);
            Candle(0,1.5,0,1,0.5);
            Candle(-4.5, 0, 0, 1, 0.5);
            Candle(4.5, 0, 0, 1, 0.5);
            break;
            
        case 6:
            rotateCircZ(0,0,0,4.5,0.2,0);
            
            break;
        case 7:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            break;
            
        case 8:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            break;
            
            
        case 9:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            rotateCircZ(0,-1,0,4.2,0.2,0);
            
            break;
            
        case 10:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            rotateCircZ(0,-1,0,4.2,0.2,0);
            
            rotateCircX(0,0,0,4.5,0.2,90);
            break;
        case 11:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            rotateCircZ(0,-1,0,4.2,0.2,0);
            
            rotateCircX(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,0,0,2,0.1,rot1);
            
            break;
            
        case 12:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            rotateCircZ(0,-1,0,4.2,0.2,0);
            
            rotateCircX(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,0,0,2,0.1,rot1);
            
            rotateCircX(0,0,0,2,0.1,rot2);
            
            
            break;
        case 13:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            rotateCircZ(0,-1,0,4.2,0.2,0);
            
            rotateCircX(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,0,0,2,0.1,rot1);
            
            rotateCircX(0,0,0,2,0.1,rot2);
            
            Ball(0,0,0, 0.5,0, 180, candle, 1);
            
            break;
            
        case 14:
            rotateCircZ(0,0,0,4.5,0.2,0);
            rotateCircZ(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,-3,0,3.2,0.2,0);
            
            rotateCircZ(0, -2, 0, 4, 0.2, 0);
            
            rotateCircZ(0,-1,0,4.2,0.2,0);
            
            rotateCircX(0,0,0,4.5,0.2,90);
            
            rotateCircZ(0,0,0,2,0.1,rot1);
            
            rotateCircX(0,0,0,2,0.1,rot2);
            
            Ball(0,0,0, 0.5,0, 180, candle, 1);
            
            chain(0,5,0);
            break;
            
        default:
            break;
    }
}

void display(){

    //  Length of axes
    const double len=25;

    //  Erase the window and the depth buffer
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    //  Set perspective
    glLoadIdentity();
    if(switchScene == 0){
        setFog();
        ox=Ex+radius*Sin(xrot);
        oz=Ez-radius*Cos(xrot);
        gluLookAt(Ex,Ey,Ez , ox,oy,oz , 0,1,0);
        //gluLookAt(Ex,Ey,Ez , 1,0,0 , 0,1,0);
        setLight(GL_LIGHT0, lx,ly,lz,ldx,ldy,ldz,sco);

        wholeMaze();
        glDisable(GL_LIGHTING);
        glDisable(GL_FOG);
    }
    else{
        double Ex1 = -2*dim*Sin(th)*Cos(ph);
        double Ey1 = +2*dim        *Sin(ph);
        double Ez1 = +2*dim*Cos(th)*Cos(ph);
        gluLookAt(Ex1,Ey1,Ez1 , 0,0,0 , 0,Cos(ph),0);
        setLight(GL_LIGHT1,10*Cos(zh), 0, 10*Sin(zh), Cos(Th)*Sin(Ph),Sin(Th)*Sin(Ph),-Cos(Ph),180);
        anotherScene();
        
        glDisable(GL_LIGHTING);
    }

    Sky(dim*6);
    //  White
   	glColor3f(1,1,1);
   	//  Draw axes
   	if (switchScene)
   	{
      	glBegin(GL_LINES);
      	glVertex3d(0.0,0.0,0.0);
      	glVertex3d(len,0.0,0.0);
      	glVertex3d(0.0,0.0,0.0);
      	glVertex3d(0.0,len,0.0);
      	glVertex3d(0.0,0.0,0.0);
      	glVertex3d(0.0,0.0,len);
      	glEnd();
      	//  Label axes
      	glRasterPos3d(len,0.0,0.0);
      	Print("X");
      	glRasterPos3d(0.0,len,0.0);
      	Print("Y");
      	glRasterPos3d(0.0,0.0,len);
      	Print("Z");
   	}

   	//  Render the scene
   	glFlush();
   	//  Make the rendered scene visible
   	glutSwapBuffers();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key,int x,int y)
{
    if(switchScene==0){
    
        if(num_points_eye == 9){
            //  Right arrow key - increase angle by 5 degrees
            if (key == GLUT_KEY_RIGHT){
                xrot += 5;

            }

            //  Left arrow key - decrease angle by 5 degrees
            else if (key == GLUT_KEY_LEFT){
                xrot -= 5;

            }
            //  Up arrow key - increase elevation by 5 degrees
            else if (key == GLUT_KEY_UP){
                Ez -= 0.2*Cos(xrot);
                Ex += 0.2*Sin(xrot);
            }

            //  Down arrow key - decrease elevation by 5 degrees
            else if (key == GLUT_KEY_DOWN){
                Ez += 0.2*Cos(xrot);
                Ex -= 0.2*Sin(xrot);
            }
        }
        Project(fov,asp,dim);
    }
    else{
        //  Right arrow key - increase angle by 5 degrees
        if (key == GLUT_KEY_RIGHT)
            th += 5;
        //  Left arrow key - decrease angle by 5 degrees
        else if (key == GLUT_KEY_LEFT)
            th -= 5;
        //  Up arrow key - increase elevation by 5 degrees
        else if (key == GLUT_KEY_UP)
            ph += 5;
        //  Down arrow key - decrease elevation by 5 degrees
        else if (key == GLUT_KEY_DOWN)
            ph -= 5;
        //  Keep angles to +/-360 degrees
        th %= 360;
        ph %= 360;
        //  Update projection
        //Project(45,asp,dim);
        Project(mode?fov:0,asp,dim);
    }



   //  Update projection
   //Project(45,asp,dim);
   
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}



/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
   //  Ratio of the width to the height of the window
   asp = (height>0) ? (double)width/height : 1;
   //  Set the viewport to the entire window
   glViewport(0,0, width,height);
   //  Set projection
   //Project(45,asp,dim);
   Project(mode?fov:0,asp,dim);
}

void eysRoutine(){
   double distance_x_e=routine[num_points_eye][0]-Ex;
   double distance_z_e=routine[num_points_eye][1]-Ez;
   //if it reaches the point
   if(-speed <= distance_x_e && distance_x_e<= speed && -speed<= distance_z_e && distance_z_e<= speed){
    Ex=routine[num_points_eye][0];
    Ez=routine[num_points_eye][1];
    num_points_eye++;
   if (num_points_eye < 8){
       double distance_x=routine[num_points_eye][0]-Ex;
       double distance_z=routine[num_points_eye][1]-Ez;
       double temp_ox=Ex+radius*Sin(xrot+TURN_RIGHT);
       double temp_oz=Ez-radius*Cos(xrot+TURN_RIGHT);
       double product=temp_ox*distance_x+temp_oz*distance_z;

       if(product > 0){
         xrot+=TURN_RIGHT;
       }
       else{
         xrot+=TURN_LEFT;
       }
    }
    else if(num_points_eye == 8){
      xrot+=TURN_LEFT;
    }
   }
   //if it doesn't reach the point.
   else{
     if(distance_x_e == 0 && distance_z_e != 0){
      if(distance_z_e > 0){
        Ez+=speed;

      }
      else{
        Ez-=speed;
      }
     }
     if(distance_x_e != 0 && distance_z_e == 0){
      if(distance_x_e > 0){
        Ex+=speed;
      }
      else{
        Ex-=speed;
      }
     }
   } 
}

void lightRoutine(){
   double distance_x_l=routine[num_points_lig][0]-lx;
   double distance_z_l=routine[num_points_lig][1]-lz;
   //if it reaches the point
   if(-speed <= distance_x_l && distance_x_l<= speed && -speed<= distance_z_l && distance_z_l<= speed){
    lx=routine[num_points_lig][0];
    lz=routine[num_points_lig][1];
    num_points_lig++;
   if (num_points_lig < 8){
       double distance_x=routine[num_points_lig][0]-lx;
       double distance_z=routine[num_points_lig][1]-lz;
       double temp_ldx=Sin(xlrot+TURN_RIGHT);
       double temp_ldz=-Cos(xlrot+TURN_RIGHT);
       double product=temp_ldx*distance_x+temp_ldz*distance_z;

       if(product > 0){
         xlrot+=TURN_RIGHT;
         ldx=temp_ldx;
         ldz=temp_ldz;
       }
       else{
         xlrot+=TURN_LEFT;
         ldx=Sin(xlrot);
         ldz=-Cos(xlrot);
       }
    }
    else if(num_points_lig == 8){
      xlrot+=TURN_LEFT;
      ldx=Sin(xlrot);
      ldz=-Cos(xlrot);
    }
   }
   //if it doesn't reach the point.
   else{
     if(distance_x_l == 0 && distance_z_l!= 0){
      if(distance_z_l > 0){
        lz+=speed;

      }
      else{
        lz-=speed;
      }
     }
     if(distance_x_l != 0 && distance_z_l == 0){
      if(distance_x_l > 0){
        lx+=speed;
      }
      else{
        lx-=speed;
      }
     }
   }   
}


/*
 *  GLUT calls this routine when the window is resized
 */
void idle()
{
   //  Elapsed time in seconds
   //double t = glutGet(GLUT_ELAPSED_TIME)/1000.0;
   
    if(switchScene == 0 && stop == 0){
        if(num_points_eye != 9){
            eysRoutine();
        }
        if(num_points_lig != 9){
            lightRoutine();
        }
    }
    //  Elapsed time in seconds
    double t = glutGet(GLUT_ELAPSED_TIME)/1000.0;
    rot1 = fmod(90*t,360.0);
    rot2= fmod(90+90*t, 360.0);
    zh = fmod(90*t,360.0);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch,int x,int y)
{
    //  Exit on ESC
    if (ch == 27)
        exit(0);

    if (ch == 't' || ch == 'T') {
        stop=1-stop;
    }
    if(ch == 'c' || ch == 'C'){
        switchScene=1-switchScene;
    }
    if(ch == 'o' || ch == 'O'){
        switchObject++;
        if(switchObject > 14){
            switchObject=0;
        }
    }
    //glutIdleFunc(stop?idle:NULL);
    
    //  Reproject
    Project(fov,asp,dim);

    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  Start up GLUT and tell it what to do
 */
int main(int argc,char* argv[])
{
	radius=sqrt((Ex-ox)*(Ex-ox)+(Ez-oz)*(Ez-oz));
   //  Initialize GLUT
   glutInit(&argc,argv);
   //  Request double buffered, true color window with Z buffering at 600x600
   glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(600,600);
   glutCreateWindow("Follow the Light!");
   //  Set callbacks
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutSpecialFunc(special);
   glutKeyboardFunc(key);
   glutIdleFunc(idle);



   //Load textures
   sky[0]=LoadTexBMP("FullMoonUp2048.bmp");
   sky[1]=LoadTexBMP("FullMoonDown2048.bmp");
   sky[2]=LoadTexBMP("FullMoonFront2048.bmp");
   sky[3]=LoadTexBMP("FullMoonRight2048.bmp");
   sky[4]=LoadTexBMP("FullMoonBack2048.bmp");
   sky[5]=LoadTexBMP("FullMoonLeft2048.bmp");
    
    bricks=LoadTexBMP("bricks.bmp");
    floors=LoadTexBMP("floor.bmp");
    metal=LoadTexBMP("gold.bmp");
    candle=LoadTexBMP("candle.bmp");
    
    //Init the displaylist for the candlesticks
    candlestick0=glGenLists (1);
    glNewList(candlestick0, GL_COMPILE);
    candlestick(-0.85,1.5,1,0.05,0.05,0.05,90);
    glEndList();
    
    candlestick1=glGenLists (1);
    glNewList(candlestick1, GL_COMPILE);
    candlestick(0.85,1.5,1,0.05,0.05,0.05,-90);
    glEndList();



   //  Pass control to GLUT so it can interact with the user
   ErrCheck("init");
   glutMainLoop();
   return 0;
}