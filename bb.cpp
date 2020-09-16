#include		<Windows.h>
#define			GL_GLEXT_PROTOTYPES
#include		<gl/gl.h>
#include		<gl/glu.h>
#include		<math.h>
#include		<vector>
#include		<sstream>
#include		<list>
#include		<queue>
#include		<map>
#include		<Shlwapi.h>
#include		<gdiplus.h>
#pragma			comment(lib, "OpenGL32.lib")
#pragma			comment(lib, "GLU32.lib")
#pragma			comment(lib, "Shlwapi.lib")
#pragma			comment(lib, "gdiplus.lib")
//#pragma			comment(lib, "Winmm.lib")//timeSetEvent, timeKillEvent

	#define		LIGHT_SYSTEM
	#define		PROFILER

int				w=0, h=0, X0=0, Y0=0;
tagPOINT		centerP, mouseP0;
tagRECT			R;
HDC__			*ghDC;
HWND__			*ghWnd;
int				kp=0;
bool			kb[256]={false};
long long		t_start=0;
unsigned		timerID=0;
std::string		path;
const unsigned char *GLversion;
int				glMajorVer, glMinorVer;
long long		broken=false;
bool			not_sent=true, noclip=true, noclipCollisionsOn=false, debuginfo=true, help=false;
unsigned long long world_seed=0;
//#ifdef _DEBUG
//int				lighting=0;
//#else
int				lighting=2;//0: fullbright, 1: blocky, 2: ambient occlusion
//#endif

#if _MSC_VER<1700//__cplusplus<201103
namespace	std
{
	inline double round(double x){return floor(x+0.5);}
	inline bool signbit(double x){return (((char*)&x)[7]>>7)!=0;}
	inline bool isnan(double x){return x!=x;}
	inline bool isinf(double x){return abs(x)==_HUGE;}
}
#endif
//Xoshiro128++	http://prng.di.unimi.it/
static unsigned seed[4];
void			xoshiro128_seed(unsigned seed){::seed[0]=seed, ::seed[1]=~seed, ::seed[2]=seed, ::seed[3]=~seed;}
void			xoshiro128_seed(unsigned v0, unsigned v1, unsigned v2, unsigned v3){seed[0]=v0, seed[1]=v1, seed[2]=v2, seed[3]=v3;}
void			xoshiro128_seed(long long lo, long long hi){*((long long*)seed)=lo, ((long long*)seed)[1]=hi;}
static inline unsigned rotl(const unsigned x, int k){return x<<k|x>>(32-k);}
unsigned		xoshiro128_next()
{
	const unsigned
		result=rotl(seed[0]+seed[3], 7)+seed[0],
		t=seed[1]<<9;

	seed[2]^=seed[0];
	seed[3]^=seed[1];
	seed[1]^=seed[2];
	seed[0]^=seed[3];

	seed[2]^=t;

	seed[3]=rotl(seed[3], 11);

	return result;
}
inline unsigned	rn24(){return xoshiro128_next()&0xFFFFFF;}

int				floor_log2(unsigned n)
{
	int logn=0;
	int sh=(((short*)&n)[1]!=0)<<4;	logn^=sh, n>>=sh;	//21.54
		sh=(((char*)&n)[1]!=0)<<3;	logn^=sh, n>>=sh;
		sh=((n&0x000000F0)!=0)<<2;	logn^=sh, n>>=sh;
		sh=((n&0x0000000C)!=0)<<1;	logn^=sh, n>>=sh;
		sh=((n&0x00000002)!=0);		logn^=sh;
	return logn;
}
inline int		maximum(int x, int y){return (x+y+abs(x-y))>>1;}
//inline float	maximum(float x, float y){return (x+y+abs(x-y))*0.5f;}
inline char		maximum(char x1, char x2, char x3, char x4, char x5, char x6)
{
	int m1=(x1+x2+abs(x1-x2)), m2=(x3+x4+abs(x3-x4)), m3=(x5+x6+abs(x5-x6));
	int m4=(m1+m2+abs(m1-m2)), m5=m3<<1;
	return (m4+m5+abs(m4-m5))>>3;
}
inline int		clamp(int lo, int x, int hi)
{
	hi<<=1;
	int temp=x+lo+abs(x-lo);
	return (temp+hi-abs(temp-hi))>>2;
}
inline float	clamp(float lo, float x, float hi)
{
	hi+=hi;
	float temp=x+lo+abs(x-lo);
	return (temp+hi-abs(temp-hi))*0.25f;
}
//inline void		clamp_vec2(float &x, float &y, float max_mag)
//{
//	float mag=sqrt(x*x+y*y), mag2=clamp(0.f, mag, max_mag)/mag;
//	x*=mag2, y*=mag2;
//}
const int		g_buf_size=1024;
int				g_buflen;
char			g_buf[g_buf_size];
void			GUIPrint(HDC__ *ghDC, int x, int y, const char *a, ...)
{
	g_buflen=vsprintf_s(g_buf, g_buf_size, a, (char*)(&a+1));
	if(g_buflen>0)
		TextOutA(ghDC, x, y, g_buf, g_buflen);
}
void			GUIPrint(HDC__ *ghDC, int x, int y, int value)
{
	g_buflen=sprintf_s(g_buf, 1024, "%d", value);
	if(g_buflen>0)
		TextOutA(ghDC, x, y, g_buf, g_buflen);
}
void			GUIVPrint(HDC__ *ghDC, int x, int y, const char *a, char *va_list)
{
	g_buflen=vsprintf_s(g_buf, g_buf_size, a, va_list);
	if(g_buflen>0)
		TextOutA(ghDC, x, y, g_buf, g_buflen);
}

void			check(int line_number)
{
	int code=glGetError();
	if(code&&!broken)
		((int*)&broken)[1]|=line_number, *(int*)&broken=code;
}
void			error(int line_number)
{
	if(!broken)
		((int*)&broken)[1]|=line_number, *(int*)&broken|=glGetError();
}

void			copy_to_clipboard(const char *a, int size)//size not including null terminator
{
	char *clipboard=(char*)LocalAlloc(LMEM_FIXED, (size+1)*sizeof(char));
	memcpy(clipboard, a, (size+1)*sizeof(char));
	clipboard[size]='\0';
	OpenClipboard(ghWnd);
	EmptyClipboard();
	SetClipboardData(CF_OEMTEXT, (void*)clipboard);
	CloseClipboard();
}
void			matrix4_to_clipboard(float *a)
{
	std::stringstream LOL_1;
	for(int k=0, idx=0;k<4;++k)
	{
		for(int k2=0;k2<4;++k2, ++idx)
			LOL_1<<a[idx]<<'\t';
		LOL_1<<"\r\n";
	}
	auto &str=LOL_1.str();
	copy_to_clipboard(str.c_str(), str.size());
}
void			bitmap_to_clipboard(int *rgb2, int w2, int h2)
{
	std::stringstream LOL_1;
	for(int ky=0;ky<h2;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			int v=0;
			switch(rgb2[w2*ky+kx]&0x00FFFFFF)
			{
			case 0x0000FF:v=1;break;
			case 0x00FF00:v=2;break;
			case 0xFF0000:v=3;break;
			}
			LOL_1<<v;
		}
		LOL_1<<"\r\n";
	}
	auto &str=LOL_1.str();
	copy_to_clipboard(str.c_str(), str.size());
}
void			vbo_to_clipboard(float *b, int size)
{
	std::stringstream LOL_1;
	for(int k=0;k<size;k+=8)
	{
		LOL_1<<k<<":\t";
		for(int k2=0;k2<8;++k2)
			LOL_1<<'\t'<<b[k+k2];
		LOL_1<<"\r\n";
	}
	auto &str=LOL_1.str();
	copy_to_clipboard(str.c_str(), str.size());
}
void			vbo_to_clipboard(std::vector<float> const &vertices, std::vector<int> const &indices, bool normals, bool texcoords)
{
	std::stringstream LOL_1;
	int stride=3+3*normals+2*texcoords;
	for(int k=0, kEnd=vertices.size();k<kEnd;k+=stride)
	{
		LOL_1<<k/stride<<":\t"<<vertices[k]<<'\t'<<vertices[k+1]<<'\t'<<vertices[k+2];
		int k2=k+3;
		if(normals)
			LOL_1<<"\t\t"<<vertices[k2]<<'\t'<<vertices[k2+1]<<'\t'<<vertices[k2+2], k2+=3;
		if(texcoords)
			LOL_1<<"\t\t"<<vertices[k2]<<'\t'<<vertices[k2+1];
		LOL_1<<"\r\n";
	}
	LOL_1<<"\r\n";
	for(int k=0, kEnd=indices.size();k<kEnd;k+=3)
		LOL_1<<k/3<<":\t"<<indices[k]<<'\t'<<indices[k+1]<<'\t'<<indices[k+2]<<"\r\n";
	auto &str=LOL_1.str();
	copy_to_clipboard(str.c_str(), str.size());
}
void			widths_to_clipboard(const char *widths)
{
	std::stringstream LOL_1;
	for(int k=0;k<128;++k)
		LOL_1<<k<<":\t"<<(int)widths[k]<<"\r\n";
	auto &str=LOL_1.str();
	copy_to_clipboard(str.c_str(), str.size());
}
template<typename T>void buffer_to_clipboard(T const *a, int size)
{
	int w=(int)sqrt((float)size), h=size/w;
	std::stringstream LOL_1;
	int idx=0;
	for(int ky=0;ky<h;++ky)
	{
		for(int kx=0;kx<w;++kx, ++idx)
			LOL_1<<a[idx]<<'\t';
		LOL_1<<"\r\n";
	}
	for(;idx<size;++idx)
		LOL_1<<a[idx]<<'\t';
	auto &str=LOL_1.str();
	copy_to_clipboard(str.c_str(), str.size());
}
typedef std::pair<std::string, double> ProfInfo;
std::vector<ProfInfo> prof;
long long		prof_t1=0;
void			prof_start()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	prof_t1=li.QuadPart;

//	prof_t1=__rdtsc();
}
void			prof_add(const char *label, int divisor=1)
{
#ifdef PROFILER
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	long long t2=li.QuadPart;
	QueryPerformanceFrequency(&li);
	prof.push_back(ProfInfo(std::string(label), 1000.*double(t2-prof_t1)/(li.QuadPart*divisor)));
	QueryPerformanceCounter(&li);
	prof_t1=li.QuadPart;

//	long long t2=__rdtsc();
//	prof.push_back(ProfInfo(std::string(label), double(t2-prof_t1)/divisor));
//	prof_t1=__rdtsc();
#endif
}
void			prof_print()
{
#ifdef PROFILER
	int xpos=w-400, xpos2=w-200;
	for(int k=0, kEnd=prof.size();k<kEnd;++k)
	{
		auto &p=prof[k];
		int ypos=k*18;
		GUIPrint(ghDC, xpos, ypos, p.first.c_str());
		GUIPrint(ghDC, xpos2, ypos, "%lf", p.second);
	//	GUIPrint(ghDC, xpos2, ypos, "%g", p.second);
	}

	//copy to clipboard
	int len=0;
	for(int k=0, kEnd=prof.size();k<kEnd;++k)
		len+=sprintf_s(g_buf+len, g_buf_size-len, "%s\t\t%lf\r\n", prof[k].first.c_str(), prof[k].second);
//	copy_to_clipboard(g_buf, len);

	//std::stringstream LOL_1;
	//for(int k=0, kEnd=prof.size();k<kEnd;++k)
	//	LOL_1<<prof[k].first<<"\t\t"<<prof[k].second<<"\r\n";
	//auto &str=LOL_1.str();
	//copy_to_clipboard(str.c_str(), str.size());

	prof.clear();
#endif
}

//vector algebra
#if 1
const float		_pi=acos(-1.f), _2pi=2*_pi, pi_2=_pi*0.5f, inv_2pi=1/_2pi, sqrt2=sqrt(2.f), torad=_pi/180, todeg=180/_pi, infinity=(float)_HUGE, inv255=1.f/255, inv256=1.f/256, inv128=1.f/128;
void			update_angle(float &angle, float &ca, float &sa){angle-=floor(angle*inv_2pi)*_2pi, ca=cos(angle), sa=sin(angle);}
void			reduce_angle(float &angle){angle-=floor(angle*inv_2pi)*_2pi;}
struct			vec2
{
	float x, y;
	vec2():x(0), y(0){}
	vec2(float x, float y):x(x), y(y){}
	void set(float x, float y){this->x=x, this->y=y;}
	vec2& operator+=(vec2 const &b){x+=b.x, y+=b.y; return *this;}
	vec2& operator-=(vec2 const &b){x-=b.x, y-=b.y; return *this;}
	vec2& operator+=(float x){this->x+=x, y+=x; return *this;}
	vec2& operator-=(float x){this->x-=x, y-=x; return *this;}
	vec2& operator*=(float x){this->x*=x, y*=x; return *this;}
	vec2& operator/=(float x){this->x/=x, y/=x; return *this;}
	float dot(vec2 const &other)const{return x*other.x+y*other.y;}
	float cross(vec2 const &other)const{return x*other.y-y*other.x;}
	float magnitude()const{return sqrt(x*x+y*y);}
	float mag_sq()const{return x*x+y*y;}
	float angle()const{return atan(y/x);}
	float angle2()const{return atan2(y, x);}
};
inline vec2		operator*(vec2 const &p, float x){return vec2(p.x*x, p.y*x);}
inline vec2		operator*(float x, vec2 const &p){return vec2(p.x*x, p.y*x);}
inline vec2		operator/(vec2 const &p, float x){return vec2(p.x/x, p.y/x);}
inline vec2		operator+(vec2 const &a, vec2 const &b){return vec2(a.x+b.x, a.y+b.y);}
inline vec2		operator-(vec2 const &a, vec2 const &b){return vec2(a.x-b.x, a.y-b.y);}
inline vec2		operator+(vec2 const &p, float x){return vec2(p.x+x, p.y+x);}
inline vec2		operator-(vec2 const &p, float x){return vec2(p.x-x, p.y-x);}
inline bool		operator==(vec2 const &a, vec2 const &b){return a.x==b.x&&a.y==b.y;}
inline bool		operator!=(vec2 const &a, vec2 const &b){return a.x!=b.x||a.y!=b.y;}
inline vec2		operator-(vec2 const &p){return vec2(-p.x, -p.y);}
struct			mat2
{
	float a, b,//row 0
		c, d;//row 1
	mat2():a(0), b(0), c(0), d(0){}
	mat2(float a, float b, float c, float d):a(a), b(b), c(c), d(d){}
};
vec2			operator*(mat2 const &m, vec2 const &v){return vec2(m.a*v.x+m.b*v.y, m.c*v.x+m.d*v.y);}
struct			vec3
{
	float x, y, z;
	vec3():x(0), y(0), z(0){}
	vec3(float x, float y, float z):x(x), y(y), z(z){}
	vec3(float gain):x(gain), y(gain), z(gain){}
	vec3(float *p):x(p[0]), y(p[1]), z(p[2]){}
	void set(float x, float y, float z){this->x=x, this->y=y, this->z=z;}
	void setzero(){x=y=z=0;}
	float& operator[](int idx){return (&x)[idx];}
	float operator[](int idx)const{return (&x)[idx];}
	vec3& operator+=(vec3 const &b){x+=b.x, y+=b.y, z+=b.z; return *this;}
	vec3& operator-=(vec3 const &b){x-=b.x, y-=b.y, z-=b.z; return *this;}
	vec3& operator+=(float x){this->x+=x, y+=x, z+=x; return *this;}
	vec3& operator-=(float x){this->x-=x, y-=x, z-=x; return *this;}
	vec3& operator*=(float x){this->x*=x, y*=x, z*=x; return *this;}
	vec3& operator/=(float x){this->x/=x, y/=x, z/=x; return *this;}
	float dot(vec3 const &other)const{return x*other.x+y*other.y+z*other.z;}
	vec3 cross(vec3 const &other)const{return vec3(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x);}
	vec3 triple_product(vec3 const &b, vec3 const &c)const;
	float magnitude()const{return sqrt(x*x+y*y+z*z);}
	float mag_sq()const{return x*x+y*y+z*z;}
	float theta()const{return atan(z/sqrt(x*x+y*y));}//vertical angle
	//float theta2()const{return atan2(y, x);}
	float phi(){return atan(y/x);}//horizontal angle
	float phi2(){return atan2(y, x);}
	bool isnan(){return x!=x||y!=y||z!=z;}
	bool isnan_or_inf(){return x!=x||y!=y||z!=z||abs(x)==infinity||abs(y)==infinity||abs(z)==infinity;}
	void clamp_xy(float max_mag)
	{
		float mag=sqrt(x*x+y*y), mag2=clamp(0.f, mag, max_mag)/mag;
		x*=mag2, y*=mag2;
	}
};
inline vec3		operator*(vec3 const &p, float x){return vec3(p.x*x, p.y*x, p.z*x);}
inline vec3		operator*(float x, vec3 const &p){return vec3(p.x*x, p.y*x, p.z*x);}
inline vec3		operator/(vec3 const &p, float x){return vec3(p.x/x, p.y/x, p.z/x);}
inline vec3		operator+(vec3 const &a, vec3 const &b){return vec3(a.x+b.x, a.y+b.y, a.z+b.z);}
inline vec3		operator-(vec3 const &a, vec3 const &b){return vec3(a.x-b.x, a.y-b.y, a.z-b.z);}
inline vec3		operator+(vec3 const &p, float x){return vec3(p.x+x, p.y+x, p.z+x);}
inline vec3		operator-(vec3 const &p, float x){return vec3(p.x-x, p.y-x, p.z-x);}
inline bool		operator==(vec3 const &a, vec3 const &b){return a.x==b.x&&a.y==b.y&&a.z==b.z;}
inline bool		operator!=(vec3 const &a, vec3 const &b){return a.x!=b.x||a.y!=b.y||a.z!=b.z;}
inline vec3		operator-(vec3 const &p){return vec3(-p.x, -p.y, -p.z);}
vec3			vec3::triple_product(vec3 const &b, vec3 const &c)const{return this->dot(c)*b-this->dot(b)*c;}
vec3			normalize(vec3 const &v){float invm=1/v.magnitude(); return invm*v;}
struct			mat3
{
	vec3	c[3];
	mat3(){}
	mat3(float gain){memset(c, 0, 9<<2), c[0][0]=c[1][1]=c[2][2]=gain;}
	mat3(vec3 const &c0, vec3 const &c1, vec3 const &c2){c[0]=c0, c[1]=c1, c[2]=c2;}
	float* data(){return &c[0][0];}
};
struct			vec4
{
	union
	{
		struct{float x, y, z, w;};
		struct{float r, g, b, a;};
	};
	vec4(){memset(&x, 0, 4<<2);}
	vec4(float gain):x(gain), y(gain), z(gain), w(gain){}
	vec4(float x, float y, float z, float w):x(x), y(y), z(z), w(w){}
	vec4(vec3 const &v, float w):x(v.x), y(v.y), z(v.z), w(w){}
	vec4(__m128 const &v){_mm_storeu_ps(&x, v);}
	float& operator[](int idx){return (&x)[idx];}
	float operator[](int idx)const{return (&x)[idx];}
	operator __m128(){return _mm_loadu_ps(&x);}
	operator vec3(){return vec3(x, y, z);}
	void set(float x, float y, float z, float w){this->x=x, this->y=y, this->z=z, this->w=w;}
	void setzero(){_mm_storeu_ps(&x, _mm_setzero_ps());}
	float dot(vec4 const &other)
	{
		__m128 a=_mm_loadu_ps(&x), b=_mm_loadu_ps(&other.x);
		a=_mm_mul_ps(a, b);
		a=_mm_hadd_ps(a, a);
		a=_mm_hadd_ps(a, a);
		float result;
		_mm_store_ss(&result, a);//does not need to be aligned
		return result;
	}
};
inline vec4		operator+(vec4 const &a, vec4 const &b)
{
	__m128 va=_mm_loadu_ps(&a.x), vb=_mm_loadu_ps(&b.x);
	va=_mm_add_ps(va, vb);
	vec4 r;
	_mm_storeu_ps(&r.x, va);
	return r;
}
inline vec4		operator*(vec4 const &v, float s)
{
	__m128 vv=_mm_loadu_ps(&v.x), vs=_mm_set1_ps(s);
	vv=_mm_mul_ps(vv, vs);
	vec4 r;
	_mm_storeu_ps(&r.x, vv);
	return r;
}
inline vec4		operator*(float s, vec4 const &v){return v*s;}
struct			mat4//column major
{
	vec4	c[4];
	mat4(){}
	mat4(const float *v){memcpy(c, v, 16<<2);}
	mat4(float gain){memset(c, 0, 16<<2), c[0][0]=c[1][1]=c[2][2]=c[3][3]=gain;}
	mat4(vec4 const &c0, vec4 const &c1, vec4 const &c2, vec4 const &c3){c[0]=c0, c[1]=c1, c[2]=c2, c[3]=c3;}
	mat4(
		float a11, float a12, float a13, float a14,
		float a21, float a22, float a23, float a24,
		float a31, float a32, float a33, float a34,
		float a41, float a42, float a43, float a44)
	{
		c[0].set(a11, a21, a31, a41);
		c[1].set(a12, a22, a32, a42);
		c[2].set(a13, a23, a33, a43);
		c[3].set(a14, a24, a34, a44);
	}
	void set(
		float a11, float a12, float a13, float a14,
		float a21, float a22, float a23, float a24,
		float a31, float a32, float a33, float a34,
		float a41, float a42, float a43, float a44)
	{
		c[0].set(a11, a21, a31, a41);
		c[1].set(a12, a22, a32, a42);
		c[2].set(a13, a23, a33, a43);
		c[3].set(a14, a24, a34, a44);
	}
	operator mat3(){return mat3((vec3)c[0], (vec3)c[1], (vec3)c[2]);}
	float* data(){return &c[0][0];}
	const float* data()const{return (float*)c;}
	void setzero(){memset(c, 0, 16<<2);}
	//void setrow(int idx, vec4 const &v){float *p=&c[0][0]+idx; p[0]=v.x, p[4]=v.y, p[8]=v.z, p[12]=v.w;}
	vec4& operator[](int idx){return c[idx];}
	vec4 operator[](int idx)const{return c[idx];}
};
inline mat4		transpose(mat4 const &m)
{
	__m128 r0=_mm_loadu_ps(&m[0][0]), r1=_mm_loadu_ps(&m[1][0]), r2=_mm_loadu_ps(&m[2][0]), r3=_mm_loadu_ps(&m[3][0]);
	_MM_TRANSPOSE4_PS(r0, r1, r2, r3);
	mat4 m2;
	_mm_storeu_ps(&m2[0][0], r0);
	_mm_storeu_ps(&m2[1][0], r1);
	_mm_storeu_ps(&m2[2][0], r2);
	_mm_storeu_ps(&m2[3][0], r3);
	return m2;
}
inline vec4		operator*(mat4 const &m, vec4 const &v)
{
	__m128 c0=_mm_loadu_ps(&m[0][0]), c1=_mm_loadu_ps(&m[1][0]), c2=_mm_loadu_ps(&m[2][0]), c3=_mm_loadu_ps(&m[3][0]);
	//	vv=_mm_loadu_ps(&v.x);
	c0=_mm_mul_ps(c0, _mm_set1_ps(v.x));
	c1=_mm_mul_ps(c1, _mm_set1_ps(v.y));
	c2=_mm_mul_ps(c2, _mm_set1_ps(v.z));
	c3=_mm_mul_ps(c3, _mm_set1_ps(v.w));
	c0=_mm_add_ps(c0, c1);
	c0=_mm_add_ps(c0, c2);
	c0=_mm_add_ps(c0, c3);
	vec4 r;
	_mm_storeu_ps(&r[0], c0);
	return r;
}
inline mat4		operator*(mat4 const &a, mat4 const &b)
{
	mat4 c=transpose(a);
	float d[]=
	{
		c[0].dot(b[0]), c[0].dot(b[1]), c[0].dot(b[2]), c[0].dot(b[3]),
		c[1].dot(b[0]), c[1].dot(b[1]), c[1].dot(b[2]), c[1].dot(b[3]),
		c[2].dot(b[0]), c[2].dot(b[1]), c[2].dot(b[2]), c[2].dot(b[3]),
		c[3].dot(b[0]), c[3].dot(b[1]), c[3].dot(b[2]), c[3].dot(b[3])
	};
	return mat4(
		vec4(c[0].dot(b[0]), c[1].dot(b[0]), c[2].dot(b[0]), c[3].dot(b[0])),
		vec4(c[0].dot(b[1]), c[1].dot(b[1]), c[2].dot(b[1]), c[3].dot(b[1])),
		vec4(c[0].dot(b[2]), c[1].dot(b[2]), c[2].dot(b[2]), c[3].dot(b[2])),
		vec4(c[0].dot(b[3]), c[1].dot(b[3]), c[2].dot(b[3]), c[3].dot(b[3])));
}
inline mat4		translate(mat4 const &m, vec3 const &delta)//from glm
{
	vec4 v2(delta, 1);
	mat4 r=m;
	r[3]=m[0]*v2[0]+m[1]*v2[1]+m[2]*v2[2]+m[3];
	return r;
}
inline mat4		rotate(mat4 const &m, float angle, vec3 const &dir)//from glm
{
	float ca=cos(angle), sa=sin(angle);
	vec3 axis=normalize(dir), temp=(1-ca)*axis;
	mat4 rotate(
		vec4(ca+temp[0]*axis[0],				temp[0]*axis[1]+sa*axis[2],		temp[0]*axis[2]-sa*axis[1],	0),//col 1
		vec4(	temp[1]*axis[0]-sa*axis[2],	ca+	temp[1]*axis[1],				temp[1]*axis[2]+sa*axis[0],	0),
		vec4(	temp[2]*axis[0]+sa*axis[1],		temp[2]*axis[1]-sa*axis[0],	ca+	temp[2]*axis[2],			0),
		vec4(0, 0, 0, 1));//col 4
	return m*rotate;
}
inline mat4		scale(mat4 const &m, vec3 const &ammount)
{
	mat4 r(m[0]*ammount.x, m[1]*ammount.y, m[2]*ammount.z, m[3]);
	return r;
}
inline mat4		lookAt(vec3 const &cam, vec3 const &center, vec3 const &up)//from glm
{
	vec3 f=normalize(center-cam),
		u=normalize(up), s=normalize(f.cross(u));
	u=s.cross(f);
	mat4 r(
		vec4(s, -s.dot(cam)),
		vec4(u, -u.dot(cam)),
		vec4(-f, f.dot(cam)),
		vec4(0, 0, 0, 1));
	r=transpose(r);
	return r;
}
inline mat4		matrixFPSViewRH(vec3 const &_cam, float pitch, float yaw)
{//https://www.3dgep.com/understanding-the-view-matrix/
	vec3 cam(_cam.y, _cam.z, _cam.x);//why yzx?
	float cos_p=cos(pitch), sin_p=sin(pitch), cos_y=cos(yaw), sin_y=sin(yaw);
	vec3
			xaxis(cos_y, 0, -sin_y),
			yaxis(sin_y*sin_p, cos_p, cos_y*sin_p),
			zaxis(sin_y*cos_p, -sin_p, cos_p*cos_y);
	
	return mat4(
		vec4(xaxis.z, yaxis.z, zaxis.z, 0),//why zxy?
		vec4(xaxis.x, yaxis.x, zaxis.x, 0),
		vec4(xaxis.y, yaxis.y, zaxis.y, 0),
		vec4(-xaxis.dot(cam), -yaxis.dot(cam), -zaxis.dot(cam), 1));
}
inline mat4		perspective(float tanfov, float ar, float znear, float zfar)
{
	return mat4(
		vec4(1/tanfov, 0, 0, 0),
		vec4(0, ar/tanfov, 0, 0),
		vec4(0, 0, -(zfar+znear)/(zfar-znear), -1),
		vec4(0, 0, -2*zfar*znear/(zfar-znear), 0));
}
struct			ivec4
{
	union
	{
		struct{int x, y, z, w;};
		struct{int x1, y1, dx, dy;};
	};
	ivec4(){_mm_storeu_si128((__m128i*)&x, _mm_setzero_si128());}
	ivec4(int x, int y, int z, int w):x(x), y(y), z(z), w(w){}
	void set(int x, int y, int z, int w){this->x=x, this->y=y, this->z=z, this->w=w;}
};
inline mat4		GetTransformInverseNoScale(const mat4 &inM)// Requires this matrix to be transform matrix, NoScale version requires this matrix be of scale 1
{
#define MakeShuffleMask(x,y,z,w)           (x | (y<<2) | (z<<4) | (w<<6))

// vec(0, 1, 2, 3) -> (vec[x], vec[y], vec[z], vec[w])
#define VecSwizzleMask(vec, mask)          _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(vec), mask))
#define VecSwizzle(vec, x, y, z, w)        VecSwizzleMask(vec, MakeShuffleMask(x,y,z,w))
#define VecSwizzle1(vec, x)                VecSwizzleMask(vec, MakeShuffleMask(x,x,x,x))
// special swizzle
#define VecSwizzle_0022(vec)               _mm_moveldup_ps(vec)
#define VecSwizzle_1133(vec)               _mm_movehdup_ps(vec)

// return (vec1[x], vec1[y], vec2[z], vec2[w])
#define VecShuffle(vec1, vec2, x,y,z,w)    _mm_shuffle_ps(vec1, vec2, MakeShuffleMask(x,y,z,w))
// special shuffle
#define VecShuffle_0101(vec1, vec2)        _mm_movelh_ps(vec1, vec2)
#define VecShuffle_2323(vec1, vec2)        _mm_movehl_ps(vec2, vec1)
	mat4 r;

	// transpose 3x3, we know m03 = m13 = m23 = 0
	__m128 t0 = VecShuffle_0101(inM[0], inM[1]); // 00, 01, 10, 11
	__m128 t1 = VecShuffle_2323(inM[0], inM[1]); // 02, 03, 12, 13
	r[0] = VecShuffle(t0, inM[2], 0,2,0,3); // 00, 10, 20, 23(=0)
	r[1] = VecShuffle(t0, inM[2], 1,3,1,3); // 01, 11, 21, 23(=0)
	r[2] = VecShuffle(t1, inM[2], 0,2,2,3); // 02, 12, 22, 23(=0)

	// last line
	r[3] =                  _mm_mul_ps(r[0], VecSwizzle1(inM[3], 0));
	r[3] = _mm_add_ps(r[3], _mm_mul_ps(r[1], VecSwizzle1(inM[3], 1)));
	r[3] = _mm_add_ps(r[3], _mm_mul_ps(r[2], VecSwizzle1(inM[3], 2)));
	r[3] = _mm_sub_ps(_mm_setr_ps(0.f, 0.f, 0.f, 1.f), r[3]);

	return r;
#undef MakeShuffleMask
#undef VecSwizzleMask
#undef VecSwizzle
#undef VecSwizzle1
#undef VecSwizzle_0022
#undef VecSwizzle_1133
#undef VecShuffle
#undef VecShuffle_0101
#undef VecShuffle_2323
}
inline mat3		normalMatrix(mat4 const &m)//inverse transpose of top left 3x3 submatrix
{
	mat4 r=GetTransformInverseNoScale(m);
	return (mat3)transpose(r);
}
#endif

float wsp=5, rsp=25;//max walking/running speed
struct			Camera
{
	const vec3	p0;
	const vec2	a0;
	const float tanfov0, dcam0, da0, mouse_sensitivity;
	vec3		p;//position
	vec2		a;//ax: yaw, ay: pitch
	float		tanfov, dcam, da, da_tfov, cax, sax, cay, say;
	Camera():p0(5, 5, 138), a0(225*torad, 324.7356103172454f*torad), tanfov0(1), dcam0(0.04f), da0(2*torad), mouse_sensitivity(0.003f),
		p(p0), a(a0), dcam(dcam0), tanfov(tanfov0), da(da0), da_tfov(tanfov), cax(cos(a.x)), sax(sin(a.x)), cay(cos(a.y)), say(sin(a.y)){}
	Camera(float x, float y, float z, float ax, float ay, float tanfov):
		p0(x, y, z), a0(ax, ay), tanfov0(tanfov), dcam0(0.04f), da0(2*torad), mouse_sensitivity(0.003f),
		p(x, y, z), a(ax, ay), tanfov(tanfov), dcam(dcam0), da(da0), da_tfov(tanfov), cax(cos(a.x)), sax(sin(a.x)), cay(cos(a.y)), say(sin(a.y)){}
	void playerRunForward	(vec3 &v, float delta){v.x=v.x+10*delta*cax, v.y+=10*delta*sax, v.clamp_xy(rsp);}
	void playerRunBack		(vec3 &v, float delta){v.x=v.x-10*delta*cax, v.y-=10*delta*sax, v.clamp_xy(rsp);}
	void playerRunLeft		(vec3 &v, float delta){v.x=v.x-10*delta*sax, v.y+=10*delta*cax, v.clamp_xy(rsp);}
	void playerRunRight		(vec3 &v, float delta){v.x=v.x+10*delta*sax, v.y-=10*delta*cax, v.clamp_xy(rsp);}
	void playerWalkForward	(vec3 &v, float delta){v.x+=delta*cax, v.y=v.y+delta*sax, v.clamp_xy(wsp);}
	void playerWalkBack		(vec3 &v, float delta){v.x-=delta*cax, v.y=v.y-delta*sax, v.clamp_xy(wsp);}
	void playerWalkLeft		(vec3 &v, float delta){v.x-=delta*sax, v.y=v.y+delta*cax, v.clamp_xy(wsp);}
	void playerWalkRight	(vec3 &v, float delta){v.x+=delta*sax, v.y=v.y-delta*cax, v.clamp_xy(wsp);}

	void moveFastForward(){p.x+=10*dcam*cax*cay,	p.y+=10*dcam*sax*cay,	p.z+=10*dcam*say;}
	void moveFastBack	(){p.x-=10*dcam*cax*cay,	p.y-=10*dcam*sax*cay,	p.z-=10*dcam*say;}
	void moveFastLeft	(){p.x-=10*dcam*sax,		p.y+=10*dcam*cax;}
	void moveFastRight	(){p.x+=10*dcam*sax,		p.y-=10*dcam*cax;}
	void moveFastUp		(){p.z+=10*dcam;}
	void moveFastDown	(){p.z-=10*dcam;}
	void moveForward	(){p.x+=dcam*cax*cay,	p.y+=dcam*sax*cay,	p.z+=dcam*say;}
	void moveBack		(){p.x-=dcam*cax*cay,	p.y-=dcam*sax*cay,	p.z-=dcam*say;}
	void moveLeft		(){p.x-=dcam*sax,		p.y+=dcam*cax;}
	void moveRight		(){p.x+=dcam*sax,		p.y-=dcam*cax;}
	void moveUp			(){p.z+=dcam;}
	void moveDown		(){p.z-=dcam;}
	void turnUp			(){update_angle(a.y+=da_tfov*da, cay, say);}
	void turnDown		(){update_angle(a.y-=da_tfov*da, cay, say);}
	void turnLeft		(){update_angle(a.x+=da_tfov*da, cax, sax);}
	void turnRight		(){update_angle(a.x-=da_tfov*da, cax, sax);}
	void turnMouse		(long lParam)
	{
		short mx=(short&)lParam, my=((short*)&lParam)[1];
		a.x+=mouse_sensitivity*da_tfov*(X0-mx), update_angle(a.x, cax, sax);
		a.y+=mouse_sensitivity*da_tfov*(Y0-my), update_angle(a.y, cay, say);
	}
	void zoomIn			(){tanfov/=1.1f, da_tfov=tanfov>1?1:tanfov;}
	void zoomOut		(){tanfov*=1.1f, da_tfov=tanfov>1?1:tanfov;}
	void reset(){p=p0, a=a0, tanfov=tanfov0,	dcam=dcam0, da_tfov=tanfov;}
	void teleport(vec3 const &p, float ax, float ay, float tanfov)
	{
		this->p=p;
		update_angle(a.x=ax, cax, sax);
		update_angle(a.y=ay, cay, say);
		this->tanfov=tanfov;
		da_tfov=tanfov>1?1:tanfov;
	}
	void faster(float A){dcam*=A;}
	void faster(){dcam*=2;}
	void slower(){dcam*=0.5;}

	vec3	viewray(){return vec3(cax*cay, sax*cay, say);}

	//vertex conversions - single precision
	void relworld2camscreen(vec3 const &d, vec3 &cp, vec2 &s)const
	{
		float cpt=d.x*cax+d.y*sax;
		cp.x=d.x*sax-d.y*cax, cp.y=cpt*say-d.z*cay, cp.z=cpt*cay+d.z*say;
		cpt=X0/(cp.z*tanfov), s.x=X0+cp.x*cpt, s.y=Y0+cp.y*cpt;
	}
	void world2camscreen(vec3 const &p_world, vec3 &cp, vec2 &s)const{relworld2camscreen(p_world-p, cp, s);}
	void relworld2cam(vec3 const &d, vec3 &cp)const
	{
		float temp=d.x*cax+d.y*sax;
		cp.x=d.x*sax-d.y*cax, cp.y=temp*say-d.z*cay, cp.z=temp*cay+d.z*say;
	}
	void world2cam(vec3 const &p, vec3 &cp)const{relworld2cam(p-this->p, cp);}
	void cam2screen(vec3 const &cp, vec2 &s)const
	{
		float temp=X0/(cp.z*tanfov);
		s.set(X0+cp.x*temp, Y0+cp.y*temp);
	}
	mat4 viewMatrix(float crouch_offset)
	{
		vec3 pos=p;
		pos.z-=crouch_offset;
		return matrixFPSViewRH(pos, a.y, a.x-_pi);
	}
} cam(5, 5, 138, 225*torad, 324.7356103172454f*torad, 1);

bool			drag=false, bypass=false, timer=false, paused=true, toSolve=true;
vec3			player_d(0.4f, 0.4f, 1.7f), player_v;
float			player_eye_height=1.6f;

enum			EngineType{E_LinearSW, E_LinearGL};
struct			Engine
{
	static int mode;
	virtual int initiate()=0;
	virtual void resize()=0;
	virtual void changeFov()=0;
	virtual void clear_screen()=0;
	//virtual void print(int x, int y, const char *a, ...)=0;
	//virtual void print(vec3 const &p, const char *a, ...)=0;
	//virtual void draw_point(vec3 const &p, int color)=0;
	//virtual void draw_line(float x1, float y1, float x2, float y2, int lineColor)=0;
	//virtual void draw_line(vec3 const &p1, vec3 const &p2, int lineColor)=0;
	//virtual void draw_triangle(vec3 const &p1, vec3 const &p2, vec3 const &p3, int color, bool transparent)=0;
	//virtual void draw_ground()=0;
	//virtual void render()=0;
	virtual void show()=0;
	//virtual void pushTexture(int *texture_)=0;
	//virtual int* popTexture()=0;
	//virtual void enqueueTextures()=0;
	//virtual void clearTextures()=0;
	virtual void finish()=0;
} *e;
int				Engine::mode=E_LinearGL;
struct			LinearSW:private Engine
{
	int			*rgb, rgbn;
	HDC__		*ghMemDC;
	HBITMAP__	*hBitmap;
	int			initiate()
	{
		rgbn=h*w;
		ghMemDC=CreateCompatibleDC(ghDC);
		tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
		hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
		SetBkMode(ghMemDC, TRANSPARENT);
		return true;
	}
	void		resize()
	{
		rgbn=h*w;
		DeleteObject((HBITMAP__*)SelectObject(ghMemDC, hBitmap));
		tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
		hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
		SetBkMode(ghMemDC, TRANSPARENT);
	}
	void		changeFov(){}
	void		clear_screen(){memset(rgb, 0xFF, rgbn*sizeof(int));}
	void		show(){BitBlt(ghDC, 0, 0, w, h, ghMemDC, 0, 0, SRCCOPY);}
	void		finish()
	{
		DeleteObject(SelectObject(ghMemDC, hBitmap));
		DeleteDC(ghMemDC);
	}
} lsw;

//OpenGL API
#if 1
bool			API_not_loaded=true;
#define			GL_MAJOR_VERSION		0x821B
#define			GL_MINOR_VERSION		0x821C
#define			GL_TEXTURE0				0x84C0
#define			GL_TEXTURE1				0x84C1
#define			GL_TEXTURE2				0x84C2
#define			GL_TEXTURE3				0x84C3
#define			GL_TEXTURE4				0x84C4
#define			GL_TEXTURE5				0x84C5
#define			GL_TEXTURE6				0x84C6
#define			GL_TEXTURE7				0x84C7
#define			GL_TEXTURE8				0x84C8
#define			GL_TEXTURE9				0x84C9
#define			GL_TEXTURE10			0x84CA
#define			GL_TEXTURE11			0x84CB
#define			GL_TEXTURE12			0x84CC
#define			GL_TEXTURE13			0x84CD
#define			GL_TEXTURE14			0x84CE
#define			GL_TEXTURE15			0x84CF
#define			GL_TEXTURE_RECTANGLE	0x84F5
#define			GL_PROGRAM_POINT_SIZE	0x8642
#define			GL_BUFFER_SIZE			0x8764
#define			GL_ARRAY_BUFFER			0x8892
#define			GL_ELEMENT_ARRAY_BUFFER	0x8893
#define			GL_STATIC_DRAW			0x88E4
#define			GL_FRAGMENT_SHADER		0x8B30
#define			GL_VERTEX_SHADER		0x8B31
#define			GL_COMPILE_STATUS		0x8B81
#define			GL_LINK_STATUS			0x8B82
#define			GL_INFO_LOG_LENGTH		0x8B84
#define			GL_DEBUG_OUTPUT			0x92E0//OpenGL 4.3+
//void			(__stdcall *glGenVertexArrays)(int n, unsigned *arrays)=nullptr;//OpenGL 3.0
//void			(__stdcall *glDeleteVertexArrays)(int n, unsigned *arrays)=nullptr;//OpenGL 3.0
void			(__stdcall *glBindVertexArray)(unsigned array)=nullptr;
void			(__stdcall *glGenBuffers)(int n, unsigned *buffers)=nullptr;
void			(__stdcall *glBindBuffer)(unsigned target, unsigned buffer)=nullptr;
void			(__stdcall *glBufferData)(unsigned target, int size, const void *data, unsigned usage)=nullptr;
void			(__stdcall *glBufferSubData)(unsigned target, int offset, int size, const void *data)=nullptr;
void			(__stdcall *glEnableVertexAttribArray)(unsigned index)=nullptr;
void			(__stdcall *glVertexAttribPointer)(unsigned index, int size, unsigned type, unsigned char normalized, int stride, const void *pointer)=nullptr;
void			(__stdcall *glDisableVertexAttribArray)(unsigned index)=nullptr;
unsigned		(__stdcall *glCreateShader)(unsigned shaderType)=nullptr;
void			(__stdcall *glShaderSource)(unsigned shader, int count, const char **string, const int *length)=nullptr;
void			(__stdcall *glCompileShader)(unsigned shader)=nullptr;
void			(__stdcall *glGetShaderiv)(unsigned shader, unsigned pname, int *params)=nullptr;
void			(__stdcall *glGetShaderInfoLog)(unsigned shader, int maxLength, int *length, char *infoLog)=nullptr;
unsigned		(__stdcall *glCreateProgram)()=nullptr;
void			(__stdcall *glAttachShader)(unsigned program, unsigned shader)=nullptr;
void			(__stdcall *glLinkProgram)(unsigned program)=nullptr;
void			(__stdcall *glGetProgramiv)(unsigned program, unsigned pname, int *params)=nullptr;
void			(__stdcall *glGetProgramInfoLog)(unsigned program, int maxLength, int *length, char *infoLog)=nullptr;
void			(__stdcall *glDetachShader)(unsigned program, unsigned shader)=nullptr;
void			(__stdcall *glDeleteShader)(unsigned shader)=nullptr;
void			(__stdcall *glUseProgram)(unsigned program)=nullptr;
int				(__stdcall *glGetAttribLocation)(unsigned program, const char *name)=nullptr;
void			(__stdcall *glDeleteProgram)(unsigned program)=nullptr;
void			(__stdcall *glDeleteBuffers)(int n, const unsigned *buffers)=nullptr;
int				(__stdcall *glGetUniformLocation)(unsigned program, const char *name)=nullptr;
void			(__stdcall *glUniform1f)(int location, float v0)=nullptr;
void			(__stdcall *glUniformMatrix3fv)(int location, int count, unsigned char transpose, const float *value)=nullptr;
void			(__stdcall *glUniformMatrix4fv)(int location, int count, unsigned char transpose, const float *value)=nullptr;
void			(__stdcall *glGetBufferParameteriv)(unsigned target, unsigned value, int *data)=nullptr;
void			(__stdcall *glActiveTexture)(unsigned texture)=nullptr;
void			(__stdcall *glUniform1i)(int location, int v0)=nullptr;
void			(__stdcall *glUniform2f)(int location, float v0, float v1)=nullptr;
void			(__stdcall *glUniform3f)(int location, float v0, float v1, float v2)=nullptr;
void			(__stdcall *glUniform3fv)(int location, int count, const float *value)=nullptr;
void			(__stdcall *glUniform4f)(int location, float v0, float v1, float v2, float v3)=nullptr;
#endif
struct			ShaderVar
{
	int *pvar;
	const char *name;
	int lineNo;//__LINE__
};
unsigned		CompileShader(const char *src, unsigned type)
{
	unsigned shaderID=glCreateShader(type);
	glShaderSource(shaderID, 1, &src, 0);
	glCompileShader(shaderID);
	int success=0;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		int infoLogLength;
		glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::vector<char> errorMessage(infoLogLength+1);
		glGetShaderInfoLog(shaderID, infoLogLength, 0, &errorMessage[0]);
		copy_to_clipboard(&errorMessage[0], infoLogLength);
		error(__LINE__);
		return 0;
	}
	return shaderID;
}
unsigned		LoadShaders(const char *vertSrc, const char *fragSrc, ShaderVar *attributes, int n_attrib, ShaderVar *uniforms, int n_unif)
{
	unsigned
		vertShaderID=CompileShader(vertSrc, GL_VERTEX_SHADER),
		fragShaderID=CompileShader(fragSrc, GL_FRAGMENT_SHADER);
//	prof_add("compile sh");
	if(!vertShaderID||!fragShaderID)
	{
		error(__LINE__);
		return 0;
	}
	unsigned ProgramID=glCreateProgram();
	glAttachShader(ProgramID, vertShaderID);
	glAttachShader(ProgramID, fragShaderID);
	glLinkProgram(ProgramID);
//	prof_add("link");
	int success=0;
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &success);
	if(!success)
	{
		int infoLogLength;
		glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &infoLogLength);
		std::vector<char> errorMessage(infoLogLength+1);
		glGetProgramInfoLog(ProgramID, infoLogLength, 0, &errorMessage[0]);
		copy_to_clipboard(&errorMessage[0], infoLogLength);
		error(__LINE__);
		return 0;
	}
	glDetachShader(ProgramID, vertShaderID);
	glDetachShader(ProgramID, fragShaderID);
	glDeleteShader(vertShaderID);
	glDeleteShader(fragShaderID);
//	prof_add("delete");
	check(__LINE__);
	for(int ka=0;ka<n_attrib;++ka)
		if((*attributes[ka].pvar=glGetAttribLocation(ProgramID, attributes[ka].name))==-1)
			error(attributes[ka].lineNo);
	for(int ku=0;ku<n_unif;++ku)
		if((*uniforms[ku].pvar=glGetUniformLocation(ProgramID, uniforms[ku].name))==-1)
			error(uniforms[ku].lineNo);
	if(broken)//
		return 0;//
	return ProgramID;
}
unsigned		current_program=0;
inline void		gl_setProgram(unsigned program)
{
	if(current_program!=program)
		glUseProgram(current_program=program);
}
void			send_color(unsigned location, int color)
{
	auto p=(unsigned char*)&color;
	static const __m128 m_255=_mm_set1_ps(inv255);
	__m128 c=_mm_castsi128_ps(_mm_set_epi32(p[3], p[2], p[1], p[0]));
	c=_mm_cvtepi32_ps(_mm_castps_si128(c));
	c=_mm_mul_ps(c, m_255);
	glUniform4f(location, c.m128_f32[0], c.m128_f32[1], c.m128_f32[2], c.m128_f32[3]);
	//glUniform4f(location, p[0]*inv255, p[1]*inv255, p[2]*inv255, p[3]*inv255);
}
void			send_color_rgb(unsigned location, int color)
{
	auto p=(unsigned char*)&color;
	static const __m128 m_255=_mm_set1_ps(inv255);
	__m128 c=_mm_castsi128_ps(_mm_set_epi32(p[3], p[2], p[1], p[0]));
	c=_mm_cvtepi32_ps(_mm_castps_si128(c));
	c=_mm_mul_ps(c, m_255);
	glUniform3f(location, c.m128_f32[0], c.m128_f32[1], c.m128_f32[2]);
}
void			send_vec3(unsigned location, vec3 const &v){glUniform3fv(location, 1, &v.x);}
void			LoadTexture(int *rgb2, int w2, int h2, unsigned &tx_id, bool alpha)
{
	glGenTextures(1, &tx_id);
	glBindTexture(GL_TEXTURE_2D, tx_id);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, alpha?GL_RGBA:GL_RGB, w2, h2, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb2);
	check(__LINE__);
}
void			select_texture(unsigned tx_id, int u_location)
{
	glActiveTexture(GL_TEXTURE0);			check(__LINE__);
	glBindTexture(GL_TEXTURE_2D, tx_id);	check(__LINE__);//select texture
	glUniform1i(u_location, 0);				check(__LINE__);
}

//BRICK BLASTER
enum			BlockType
{
	B_AIR, B_GLASS,//air is zero, last transparent block is glass
	B_TORCH,
	B_DIRT, B_WALL, B_BEDROCK,
		ALL_BLOCKS
};
const int		nTextures=4, nTransparentTextures=1;
unsigned		textures[nTextures+1]={0};
void			generate_textures()
{
	int txw=16, txh=16, tx_size=txw*txh, *rgb=(int*)malloc(tx_size<<2);//0xAABBGGRR

	for(int ky=0;ky<txh;++ky)//bedrock - grey lines
	{
		int color;
		auto p=(unsigned char*)&color;
		p[0]=p[1]=p[2]=rand();
		p[3]=0xFF;
		std::fill(rgb+txw*ky, rgb+txw*(ky+1), color);
	}
	LoadTexture(rgb, txw, txh, textures[B_BEDROCK], true);

	for(int k=0;k<tx_size;++k)//dirt - brown
	{
		auto p=(unsigned char*)(rgb+k);
		p[0]=0x7F+(rand()&0x1F);
	//	p[0]=rand()&0x7F;
	//	p[0]=0x7F;
		p[1]=p[0]*2/3;
	//	p[1]=(rand()&0xFF)>>1;
		p[2]=p[0]/3;
		p[3]=0xFF;
	}
	LoadTexture(rgb, txw, txh, textures[B_DIRT], true);

	for(int k=0;k<tx_size;++k)//wall - grey
	{
		auto p=(unsigned char*)(rgb+k);
		p[0]=p[1]=p[2]=rand();
		p[3]=0xFF;
	}
	LoadTexture(rgb, txw, txh, textures[B_WALL], true);

	for(int k=0;k<tx_size;++k)//torch - red
	{
		auto p=(unsigned char*)(rgb+k);
		p[0]=rand();
		p[1]=p[2]=0;
		p[3]=0xFF;
	}
	LoadTexture(rgb, txw, txh, textures[B_TORCH], true);

	for(int k=0;k<tx_size;++k)//glass - small or zero alpha
	{
		rgb[k]=rand()<<15|rand();
		auto p=(unsigned char*)(rgb+k);
		p[3]<<=2;
	//	p[3]>>=2;
	//	p[3]=0xFF;
	}
	LoadTexture(rgb, txw, txh, textures[B_GLASS], true);

	free(rgb);
}
//OpenSimplex Noise		2014-09-19		99.77ms
#if 1//https://gist.github.com/KdotJPG/b1270127455a94ac5d19
const double STRETCH_CONSTANT_2D=-0.211324865405187;//(1/sqrt(2+1)-1)/2;
const double SQUISH_CONSTANT_2D=0.366025403784439;	//(sqrt(2+1)-1)/2;
const double STRETCH_CONSTANT_3D=-1./6;				//(1/sqrt(3+1)-1)/3;
const double SQUISH_CONSTANT_3D=1./3;				//(sqrt(3+1)-1)/3;
//const double STRETCH_CONSTANT_4D=-0.138196601125011;//(1/sqrt(4+1)-1)/4;
//const double SQUISH_CONSTANT_4D=0.309016994374947;	//(sqrt(4+1)-1)/4;

const double NORM_CONSTANT_2D=47;
const double NORM_CONSTANT_3D=103;
//const double NORM_CONSTANT_4D=30;
	
const long long DEFAULT_SEED=0;
class			OpenSimplexNoise
{
	static const char gradients3D[72], gradients2D[16];
	short perm[256], permGradIndex3D[256];
	double extrapolate(int xsb, int ysb, double dx, double dy)
	{
		int index=perm[(perm[xsb&0xFF]+ysb)&0xFF]&0x0E;
		return gradients2D[index]*dx+gradients2D[index+1]*dy;
	}
	double extrapolate(int xsb, int ysb, int zsb, double dx, double dy, double dz)
	{
		int index=permGradIndex3D[(perm[(perm[xsb&0xFF]+ysb)&0xFF]+zsb)&0xFF];
		return gradients3D[index]*dx+ gradients3D[index+1]*dy+gradients3D[index+2]*dz;
	}
public:
	//Initializes the class using a permutation array generated from a 64-bit seed.
	//Generates a proper permutation (i.e. doesn't merely perform N successive pair swaps on a base array)
	//Uses a simple 64-bit LCG.
	void seed(long long seed)
	{
	//	perm = new short[256];
	//	permGradIndex3D = new short[256];
		short source[256];
		for(short i=0;i<256;++i)
			source[i]=i;
		seed=seed*6364136223846793005+1442695040888963407;
		seed=seed*6364136223846793005+1442695040888963407;
		seed=seed*6364136223846793005+1442695040888963407;
		for(int i=255;i>=0;--i)
		{
			seed=seed*6364136223846793005l+1442695040888963407l;
			int r=(int)((seed+31)%(i+1));
			if(r<0)
				r+=i+1;
			perm[i]=source[r];
			permGradIndex3D[i]=(short)((perm[i]%(sizeof(gradients3D)/3))*3);
			source[r] = source[i];
		}
	}
	double evaluate(double x, double y)//2D OpenSimplex Noise.
	{
		double stretchOffset=(x+y)*STRETCH_CONSTANT_2D;//Place input coordinates onto grid.
		double xs=x+stretchOffset, ys=y+stretchOffset;
		
		int xsb=(int)floor(xs), ysb=(int)floor(ys);//Floor to get grid coordinates of rhombus (stretched square) super-cell origin.
		
		double squishOffset=(xsb+ysb)*SQUISH_CONSTANT_2D;//Skew out to get actual coordinates of rhombus origin. We'll need these later.
		double xb=xsb+squishOffset, yb=ysb+squishOffset;
		
		double xins=xs-xsb, yins=ys-ysb;//Compute grid coordinates relative to rhombus origin.
		double inSum=xins+yins;//Sum those together to get a value that determines which region we're in.
		double dx0=x-xb, dy0=y-yb;//Positions relative to origin point.
		
		double dx_ext, dy_ext;//We'll be defining these inside the next block and using them afterwards.
		int xsv_ext, ysv_ext;
		double value = 0;

		double dx1=dx0-1-SQUISH_CONSTANT_2D;//Contribution (1,0)
		double dy1=dy0-0-SQUISH_CONSTANT_2D;
		double attn1=2-dx1*dx1-dy1*dy1;
		if(attn1>0)
		{
			attn1*=attn1;
			value+=attn1*attn1*extrapolate(xsb+1, ysb+0, dx1, dy1);
		}

		double dx2=dx0-0-SQUISH_CONSTANT_2D;//Contribution (0,1)
		double dy2=dy0-1-SQUISH_CONSTANT_2D;
		double attn2=2-dx2*dx2-dy2*dy2;
		if(attn2>0)
		{
			attn2*=attn2;
			value+=attn2*attn2*extrapolate(xsb+0, ysb+1, dx2, dy2);
		}
		
		if (inSum <= 1)//We're inside the triangle (2-Simplex) at (0,0)
		{
			double zins=1-inSum;
			if(zins>xins||zins>yins)//(0,0) is one of the closest two triangular vertices
			{
				if (xins > yins)
				{
					xsv_ext=xsb+1, ysv_ext=ysb-1;
					dx_ext=dx0-1, dy_ext=dy0+1;
				}
				else
				{
					xsv_ext=xsb-1, ysv_ext=ysb+1;
					dx_ext=dx0+1, dy_ext=dy0-1;
				}
			}
			else//(1,0) and (0,1) are the closest two vertices.
			{
				xsv_ext=xsb+1, ysv_ext=ysb+1;
				dx_ext=dx0-1-2*SQUISH_CONSTANT_2D, dy_ext=dy0-1-2*SQUISH_CONSTANT_2D;
			}
		}
		else//We're inside the triangle (2-Simplex) at (1,1)
		{
			double zins=2-inSum;
			if (zins<xins||zins<yins)//(0,0) is one of the closest two triangular vertices
			{
				if(xins>yins)
				{
					xsv_ext=xsb+2, ysv_ext=ysb+0;
					dx_ext=dx0-2-2*SQUISH_CONSTANT_2D, dy_ext=dy0+0-2*SQUISH_CONSTANT_2D;
				}
				else
				{
					xsv_ext=xsb+0;
					ysv_ext=ysb+2;
					dx_ext=dx0+0-2*SQUISH_CONSTANT_2D;
					dy_ext=dy0-2-2*SQUISH_CONSTANT_2D;
				}
			}
			else//(1,0) and (0,1) are the closest two vertices.
			{
				dx_ext=dx0, dy_ext=dy0;
				xsv_ext=xsb, ysv_ext=ysb;
			}
			++xsb, ++ysb;
			dx0=dx0-1-2*SQUISH_CONSTANT_2D, dy0=dy0-1-2*SQUISH_CONSTANT_2D;
		}
		double attn0=2-dx0*dx0-dy0*dy0;//Contribution (0,0) or (1,1)
		if(attn0>0)
		{
			attn0*=attn0;
			value+=attn0*attn0*extrapolate(xsb, ysb, dx0, dy0);
		}
		double attn_ext=2-dx_ext*dx_ext-dy_ext*dy_ext;//Extra Vertex
		if(attn_ext>0)
		{
			attn_ext*=attn_ext;
			value+=attn_ext*attn_ext*extrapolate(xsv_ext, ysv_ext, dx_ext, dy_ext);
		}
		return value/NORM_CONSTANT_2D;
	}
	double evaluate(double x, double y, double z)//3D OpenSimplex Noise.
	{
	
		//Place input coordinates on simplectic honeycomb.
		double stretchOffset = (x + y + z) * STRETCH_CONSTANT_3D;
		double xs = x + stretchOffset;
		double ys = y + stretchOffset;
		double zs = z + stretchOffset;
		
		//Floor to get simplectic honeycomb coordinates of rhombohedron (stretched cube) super-cell origin.
		int xsb = (int)floor(xs);
		int ysb = (int)floor(ys);
		int zsb = (int)floor(zs);
		
		//Skew out to get actual coordinates of rhombohedron origin. We'll need these later.
		double squishOffset = (xsb + ysb + zsb) * SQUISH_CONSTANT_3D;
		double xb = xsb + squishOffset;
		double yb = ysb + squishOffset;
		double zb = zsb + squishOffset;
		
		//Compute simplectic honeycomb coordinates relative to rhombohedral origin.
		double xins = xs - xsb;
		double yins = ys - ysb;
		double zins = zs - zsb;
		
		//Sum those together to get a value that determines which region we're in.
		double inSum = xins + yins + zins;

		//Positions relative to origin point.
		double dx0 = x - xb;
		double dy0 = y - yb;
		double dz0 = z - zb;
		
		//We'll be defining these inside the next block and using them afterwards.
		double dx_ext0, dy_ext0, dz_ext0;
		double dx_ext1, dy_ext1, dz_ext1;
		int xsv_ext0, ysv_ext0, zsv_ext0;
		int xsv_ext1, ysv_ext1, zsv_ext1;
		
		double value = 0;
		if (inSum <= 1) { //We're inside the tetrahedron (3-Simplex) at (0,0,0)
			
			//Determine which two of (0,0,1), (0,1,0), (1,0,0) are closest.
			byte aPoint = 0x01;
			double aScore = xins;
			byte bPoint = 0x02;
			double bScore = yins;
			if (aScore >= bScore && zins > bScore) {
				bScore = zins;
				bPoint = 0x04;
			} else if (aScore < bScore && zins > aScore) {
				aScore = zins;
				aPoint = 0x04;
			}
			
			//Now we determine the two lattice points not part of the tetrahedron that may contribute.
			//This depends on the closest two tetrahedral vertices, including (0,0,0)
			double wins = 1 - inSum;
			if (wins > aScore || wins > bScore) { //(0,0,0) is one of the closest two tetrahedral vertices.
				byte c = (bScore > aScore ? bPoint : aPoint); //Our other closest vertex is the closest out of a and b.
				
				if ((c & 0x01) == 0) {
					xsv_ext0 = xsb - 1;
					xsv_ext1 = xsb;
					dx_ext0 = dx0 + 1;
					dx_ext1 = dx0;
				} else {
					xsv_ext0 = xsv_ext1 = xsb + 1;
					dx_ext0 = dx_ext1 = dx0 - 1;
				}

				if ((c & 0x02) == 0) {
					ysv_ext0 = ysv_ext1 = ysb;
					dy_ext0 = dy_ext1 = dy0;
					if ((c & 0x01) == 0) {
						ysv_ext1 -= 1;
						dy_ext1 += 1;
					} else {
						ysv_ext0 -= 1;
						dy_ext0 += 1;
					}
				} else {
					ysv_ext0 = ysv_ext1 = ysb + 1;
					dy_ext0 = dy_ext1 = dy0 - 1;
				}

				if ((c & 0x04) == 0) {
					zsv_ext0 = zsb;
					zsv_ext1 = zsb - 1;
					dz_ext0 = dz0;
					dz_ext1 = dz0 + 1;
				} else {
					zsv_ext0 = zsv_ext1 = zsb + 1;
					dz_ext0 = dz_ext1 = dz0 - 1;
				}
			} else { //(0,0,0) is not one of the closest two tetrahedral vertices.
				byte c = (byte)(aPoint | bPoint); //Our two extra vertices are determined by the closest two.
				
				if ((c & 0x01) == 0) {
					xsv_ext0 = xsb;
					xsv_ext1 = xsb - 1;
					dx_ext0 = dx0 - 2 * SQUISH_CONSTANT_3D;
					dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D;
				} else {
					xsv_ext0 = xsv_ext1 = xsb + 1;
					dx_ext0 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D;
					dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D;
				}

				if ((c & 0x02) == 0) {
					ysv_ext0 = ysb;
					ysv_ext1 = ysb - 1;
					dy_ext0 = dy0 - 2 * SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D;
				} else {
					ysv_ext0 = ysv_ext1 = ysb + 1;
					dy_ext0 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D;
				}

				if ((c & 0x04) == 0) {
					zsv_ext0 = zsb;
					zsv_ext1 = zsb - 1;
					dz_ext0 = dz0 - 2 * SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D;
				} else {
					zsv_ext0 = zsv_ext1 = zsb + 1;
					dz_ext0 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D;
				}
			}

			//Contribution (0,0,0)
			double attn0 = 2 - dx0 * dx0 - dy0 * dy0 - dz0 * dz0;
			if (attn0 > 0) {
				attn0 *= attn0;
				value += attn0 * attn0 * extrapolate(xsb + 0, ysb + 0, zsb + 0, dx0, dy0, dz0);
			}

			//Contribution (1,0,0)
			double dx1 = dx0 - 1 - SQUISH_CONSTANT_3D;
			double dy1 = dy0 - 0 - SQUISH_CONSTANT_3D;
			double dz1 = dz0 - 0 - SQUISH_CONSTANT_3D;
			double attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1;
			if (attn1 > 0) {
				attn1 *= attn1;
				value += attn1 * attn1 * extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1);
			}

			//Contribution (0,1,0)
			double dx2 = dx0 - 0 - SQUISH_CONSTANT_3D;
			double dy2 = dy0 - 1 - SQUISH_CONSTANT_3D;
			double dz2 = dz1;
			double attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2;
			if (attn2 > 0) {
				attn2 *= attn2;
				value += attn2 * attn2 * extrapolate(xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2);
			}

			//Contribution (0,0,1)
			double dx3 = dx2;
			double dy3 = dy1;
			double dz3 = dz0 - 1 - SQUISH_CONSTANT_3D;
			double attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3;
			if (attn3 > 0) {
				attn3 *= attn3;
				value += attn3 * attn3 * extrapolate(xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3);
			}
		} else if (inSum >= 2) { //We're inside the tetrahedron (3-Simplex) at (1,1,1)
		
			//Determine which two tetrahedral vertices are the closest, out of (1,1,0), (1,0,1), (0,1,1) but not (1,1,1).
			byte aPoint = 0x06;
			double aScore = xins;
			byte bPoint = 0x05;
			double bScore = yins;
			if (aScore <= bScore && zins < bScore) {
				bScore = zins;
				bPoint = 0x03;
			} else if (aScore > bScore && zins < aScore) {
				aScore = zins;
				aPoint = 0x03;
			}
			
			//Now we determine the two lattice points not part of the tetrahedron that may contribute.
			//This depends on the closest two tetrahedral vertices, including (1,1,1)
			double wins = 3 - inSum;
			if (wins < aScore || wins < bScore) { //(1,1,1) is one of the closest two tetrahedral vertices.
				byte c = (bScore < aScore ? bPoint : aPoint); //Our other closest vertex is the closest out of a and b.
				
				if ((c & 0x01) != 0) {
					xsv_ext0 = xsb + 2;
					xsv_ext1 = xsb + 1;
					dx_ext0 = dx0 - 2 - 3 * SQUISH_CONSTANT_3D;
					dx_ext1 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D;
				} else {
					xsv_ext0 = xsv_ext1 = xsb;
					dx_ext0 = dx_ext1 = dx0 - 3 * SQUISH_CONSTANT_3D;
				}

				if ((c & 0x02) != 0) {
					ysv_ext0 = ysv_ext1 = ysb + 1;
					dy_ext0 = dy_ext1 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D;
					if ((c & 0x01) != 0) {
						ysv_ext1 += 1;
						dy_ext1 -= 1;
					} else {
						ysv_ext0 += 1;
						dy_ext0 -= 1;
					}
				} else {
					ysv_ext0 = ysv_ext1 = ysb;
					dy_ext0 = dy_ext1 = dy0 - 3 * SQUISH_CONSTANT_3D;
				}

				if ((c & 0x04) != 0) {
					zsv_ext0 = zsb + 1;
					zsv_ext1 = zsb + 2;
					dz_ext0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 2 - 3 * SQUISH_CONSTANT_3D;
				} else {
					zsv_ext0 = zsv_ext1 = zsb;
					dz_ext0 = dz_ext1 = dz0 - 3 * SQUISH_CONSTANT_3D;
				}
			} else { //(1,1,1) is not one of the closest two tetrahedral vertices.
				byte c = (byte)(aPoint & bPoint); //Our two extra vertices are determined by the closest two.
				
				if ((c & 0x01) != 0) {
					xsv_ext0 = xsb + 1;
					xsv_ext1 = xsb + 2;
					dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D;
					dx_ext1 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D;
				} else {
					xsv_ext0 = xsv_ext1 = xsb;
					dx_ext0 = dx0 - SQUISH_CONSTANT_3D;
					dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
				}

				if ((c & 0x02) != 0) {
					ysv_ext0 = ysb + 1;
					ysv_ext1 = ysb + 2;
					dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D;
				} else {
					ysv_ext0 = ysv_ext1 = ysb;
					dy_ext0 = dy0 - SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
				}

				if ((c & 0x04) != 0) {
					zsv_ext0 = zsb + 1;
					zsv_ext1 = zsb + 2;
					dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D;
				} else {
					zsv_ext0 = zsv_ext1 = zsb;
					dz_ext0 = dz0 - SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
				}
			}
			
			//Contribution (1,1,0)
			double dx3 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D;
			double dy3 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D;
			double dz3 = dz0 - 0 - 2 * SQUISH_CONSTANT_3D;
			double attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3;
			if (attn3 > 0) {
				attn3 *= attn3;
				value += attn3 * attn3 * extrapolate(xsb + 1, ysb + 1, zsb + 0, dx3, dy3, dz3);
			}

			//Contribution (1,0,1)
			double dx2 = dx3;
			double dy2 = dy0 - 0 - 2 * SQUISH_CONSTANT_3D;
			double dz2 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D;
			double attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2;
			if (attn2 > 0) {
				attn2 *= attn2;
				value += attn2 * attn2 * extrapolate(xsb + 1, ysb + 0, zsb + 1, dx2, dy2, dz2);
			}

			//Contribution (0,1,1)
			double dx1 = dx0 - 0 - 2 * SQUISH_CONSTANT_3D;
			double dy1 = dy3;
			double dz1 = dz2;
			double attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1;
			if (attn1 > 0) {
				attn1 *= attn1;
				value += attn1 * attn1 * extrapolate(xsb + 0, ysb + 1, zsb + 1, dx1, dy1, dz1);
			}

			//Contribution (1,1,1)
			dx0 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D;
			dy0 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D;
			dz0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D;
			double attn0 = 2 - dx0 * dx0 - dy0 * dy0 - dz0 * dz0;
			if (attn0 > 0) {
				attn0 *= attn0;
				value += attn0 * attn0 * extrapolate(xsb + 1, ysb + 1, zsb + 1, dx0, dy0, dz0);
			}
		} else { //We're inside the octahedron (Rectified 3-Simplex) in between.
			double aScore;
			byte aPoint;
			boolean aIsFurtherSide;
			double bScore;
			byte bPoint;
			boolean bIsFurtherSide;

			//Decide between point (0,0,1) and (1,1,0) as closest
			double p1 = xins + yins;
			if (p1 > 1) {
				aScore = p1 - 1;
				aPoint = 0x03;
				aIsFurtherSide = true;
			} else {
				aScore = 1 - p1;
				aPoint = 0x04;
				aIsFurtherSide = false;
			}

			//Decide between point (0,1,0) and (1,0,1) as closest
			double p2 = xins + zins;
			if (p2 > 1) {
				bScore = p2 - 1;
				bPoint = 0x05;
				bIsFurtherSide = true;
			} else {
				bScore = 1 - p2;
				bPoint = 0x02;
				bIsFurtherSide = false;
			}
			
			//The closest out of the two (1,0,0) and (0,1,1) will replace the furthest out of the two decided above, if closer.
			double p3 = yins + zins;
			if (p3 > 1) {
				double score = p3 - 1;
				if (aScore <= bScore && aScore < score) {
					aScore = score;
					aPoint = 0x06;
					aIsFurtherSide = true;
				} else if (aScore > bScore && bScore < score) {
					bScore = score;
					bPoint = 0x06;
					bIsFurtherSide = true;
				}
			} else {
				double score = 1 - p3;
				if (aScore <= bScore && aScore < score) {
					aScore = score;
					aPoint = 0x01;
					aIsFurtherSide = false;
				} else if (aScore > bScore && bScore < score) {
					bScore = score;
					bPoint = 0x01;
					bIsFurtherSide = false;
				}
			}
			
			//Where each of the two closest points are determines how the extra two vertices are calculated.
			if (aIsFurtherSide == bIsFurtherSide) {
				if (aIsFurtherSide) { //Both closest points on (1,1,1) side

					//One of the two extra points is (1,1,1)
					dx_ext0 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D;
					dy_ext0 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D;
					dz_ext0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D;
					xsv_ext0 = xsb + 1;
					ysv_ext0 = ysb + 1;
					zsv_ext0 = zsb + 1;

					//Other extra point is based on the shared axis.
					byte c = (byte)(aPoint & bPoint);
					if ((c & 0x01) != 0) {
						dx_ext1 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D;
						dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
						dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
						xsv_ext1 = xsb + 2;
						ysv_ext1 = ysb;
						zsv_ext1 = zsb;
					} else if ((c & 0x02) != 0) {
						dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
						dy_ext1 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D;
						dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
						xsv_ext1 = xsb;
						ysv_ext1 = ysb + 2;
						zsv_ext1 = zsb;
					} else {
						dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
						dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
						dz_ext1 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D;
						xsv_ext1 = xsb;
						ysv_ext1 = ysb;
						zsv_ext1 = zsb + 2;
					}
				} else {//Both closest points on (0,0,0) side

					//One of the two extra points is (0,0,0)
					dx_ext0 = dx0;
					dy_ext0 = dy0;
					dz_ext0 = dz0;
					xsv_ext0 = xsb;
					ysv_ext0 = ysb;
					zsv_ext0 = zsb;

					//Other extra point is based on the omitted axis.
					byte c = (byte)(aPoint | bPoint);
					if ((c & 0x01) == 0) {
						dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D;
						dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D;
						dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D;
						xsv_ext1 = xsb - 1;
						ysv_ext1 = ysb + 1;
						zsv_ext1 = zsb + 1;
					} else if ((c & 0x02) == 0) {
						dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D;
						dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D;
						dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D;
						xsv_ext1 = xsb + 1;
						ysv_ext1 = ysb - 1;
						zsv_ext1 = zsb + 1;
					} else {
						dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D;
						dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D;
						dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D;
						xsv_ext1 = xsb + 1;
						ysv_ext1 = ysb + 1;
						zsv_ext1 = zsb - 1;
					}
				}
			} else { //One point on (0,0,0) side, one point on (1,1,1) side
				byte c1, c2;
				if (aIsFurtherSide) {
					c1 = aPoint;
					c2 = bPoint;
				} else {
					c1 = bPoint;
					c2 = aPoint;
				}

				//One contribution is a permutation of (1,1,-1)
				if ((c1 & 0x01) == 0) {
					dx_ext0 = dx0 + 1 - SQUISH_CONSTANT_3D;
					dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D;
					dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D;
					xsv_ext0 = xsb - 1;
					ysv_ext0 = ysb + 1;
					zsv_ext0 = zsb + 1;
				} else if ((c1 & 0x02) == 0) {
					dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D;
					dy_ext0 = dy0 + 1 - SQUISH_CONSTANT_3D;
					dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D;
					xsv_ext0 = xsb + 1;
					ysv_ext0 = ysb - 1;
					zsv_ext0 = zsb + 1;
				} else {
					dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D;
					dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D;
					dz_ext0 = dz0 + 1 - SQUISH_CONSTANT_3D;
					xsv_ext0 = xsb + 1;
					ysv_ext0 = ysb + 1;
					zsv_ext0 = zsb - 1;
				}

				//One contribution is a permutation of (0,0,2)
				dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
				dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
				dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
				xsv_ext1 = xsb;
				ysv_ext1 = ysb;
				zsv_ext1 = zsb;
				if ((c2 & 0x01) != 0) {
					dx_ext1 -= 2;
					xsv_ext1 += 2;
				} else if ((c2 & 0x02) != 0) {
					dy_ext1 -= 2;
					ysv_ext1 += 2;
				} else {
					dz_ext1 -= 2;
					zsv_ext1 += 2;
				}
			}

			//Contribution (1,0,0)
			double dx1 = dx0 - 1 - SQUISH_CONSTANT_3D;
			double dy1 = dy0 - 0 - SQUISH_CONSTANT_3D;
			double dz1 = dz0 - 0 - SQUISH_CONSTANT_3D;
			double attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1;
			if (attn1 > 0) {
				attn1 *= attn1;
				value += attn1 * attn1 * extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1);
			}

			//Contribution (0,1,0)
			double dx2 = dx0 - 0 - SQUISH_CONSTANT_3D;
			double dy2 = dy0 - 1 - SQUISH_CONSTANT_3D;
			double dz2 = dz1;
			double attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2;
			if (attn2 > 0) {
				attn2 *= attn2;
				value += attn2 * attn2 * extrapolate(xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2);
			}

			//Contribution (0,0,1)
			double dx3 = dx2;
			double dy3 = dy1;
			double dz3 = dz0 - 1 - SQUISH_CONSTANT_3D;
			double attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3;
			if (attn3 > 0) {
				attn3 *= attn3;
				value += attn3 * attn3 * extrapolate(xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3);
			}

			//Contribution (1,1,0)
			double dx4 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D;
			double dy4 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D;
			double dz4 = dz0 - 0 - 2 * SQUISH_CONSTANT_3D;
			double attn4 = 2 - dx4 * dx4 - dy4 * dy4 - dz4 * dz4;
			if (attn4 > 0) {
				attn4 *= attn4;
				value += attn4 * attn4 * extrapolate(xsb + 1, ysb + 1, zsb + 0, dx4, dy4, dz4);
			}

			//Contribution (1,0,1)
			double dx5 = dx4;
			double dy5 = dy0 - 0 - 2 * SQUISH_CONSTANT_3D;
			double dz5 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D;
			double attn5 = 2 - dx5 * dx5 - dy5 * dy5 - dz5 * dz5;
			if (attn5 > 0) {
				attn5 *= attn5;
				value += attn5 * attn5 * extrapolate(xsb + 1, ysb + 0, zsb + 1, dx5, dy5, dz5);
			}

			//Contribution (0,1,1)
			double dx6 = dx0 - 0 - 2 * SQUISH_CONSTANT_3D;
			double dy6 = dy4;
			double dz6 = dz5;
			double attn6 = 2 - dx6 * dx6 - dy6 * dy6 - dz6 * dz6;
			if (attn6 > 0) {
				attn6 *= attn6;
				value += attn6 * attn6 * extrapolate(xsb + 0, ysb + 1, zsb + 1, dx6, dy6, dz6);
			}
		}
 
		//First extra vertex
		double attn_ext0 = 2 - dx_ext0 * dx_ext0 - dy_ext0 * dy_ext0 - dz_ext0 * dz_ext0;
		if (attn_ext0 > 0)
		{
			attn_ext0 *= attn_ext0;
			value += attn_ext0 * attn_ext0 * extrapolate(xsv_ext0, ysv_ext0, zsv_ext0, dx_ext0, dy_ext0, dz_ext0);
		}

		//Second extra vertex
		double attn_ext1 = 2 - dx_ext1 * dx_ext1 - dy_ext1 * dy_ext1 - dz_ext1 * dz_ext1;
		if (attn_ext1 > 0)
		{
			attn_ext1 *= attn_ext1;
			value += attn_ext1 * attn_ext1 * extrapolate(xsv_ext1, ysv_ext1, zsv_ext1, dx_ext1, dy_ext1, dz_ext1);
		}
		
		return value / NORM_CONSTANT_3D;
	}
};
//Gradients for 2D. They approximate the directions to the vertices of an octagon from the center.
const char OpenSimplexNoise::gradients2D[]=
{
	 5,  2,    2,  5,
	-5,  2,   -2,  5,
	 5, -2,    2, -5,
	-5, -2,   -2, -5,
};
//Gradients for 3D. They approximate the directions to the vertices of a rhombicuboctahedron from the center, skewed so
//that the triangular and square facets can be inscribed inside circles of the same radius.
const char OpenSimplexNoise::gradients3D[]=
{
	-11,  4,  4,     -4,  11,  4,    -4,  4,  11,
	 11,  4,  4,      4,  11,  4,     4,  4,  11,
	-11, -4,  4,     -4, -11,  4,    -4, -4,  11,
	 11, -4,  4,      4, -11,  4,     4, -4,  11,
	-11,  4, -4,     -4,  11, -4,    -4,  4, -11,
	 11,  4, -4,      4,  11, -4,     4,  4, -11,
	-11, -4, -4,     -4, -11, -4,    -4, -4, -11,
	 11, -4, -4,      4, -11, -4,     4, -4, -11,
};
OpenSimplexNoise pdf_surface;
//OpenSimplexNoise pdf_caves, pdf_surface;
#endif
 //Perlin Noise		1983				24.30ms
#if 1
#define			mix(x, y, ay)	((x)+((y)-(x))*(ay))
struct			PerlinNoise2D//within +-sqrt(2)/2=+-0.707
{
	static const int log_g=3, g_side=1<<log_g, g_size=1<<(log_g<<1);
	vec2		normals[g_size];
	std::hash<int> hash;
	void generate_normals()
	{
		const float inv_rmax=2.f/0xFFFFFF;
		for(int k=0;k<g_size;++k)
		{
			auto &n=normals[k];
			do
			{
				n.x=rn24()*inv_rmax-1;
				n.y=rn24()*inv_rmax-1;
			}
			while(n.mag_sq()>1);
			n*=1.f/n.magnitude();
		}
	}
	void		seed(unsigned seed)
	{
		srand(seed);
		xoshiro128_seed(rand()<<15|rand(), rand()<<15|rand(), rand()<<15|rand(), rand()<<15|rand());
		generate_normals();
	}
	void		seed(unsigned long long seed)
	{
		srand((unsigned)seed);
		unsigned long long s0=(long long)rand()<<45|(long long)rand()<<30|rand()<<15|rand();
		srand((unsigned)(seed>>32));
		unsigned long long s1=(long long)rand()<<45|(long long)rand()<<30|rand()<<15|rand();
		xoshiro128_seed(s0, s1);
		generate_normals();
	}
	unsigned hash3(int x)
	{
		unsigned h=hash(x);
		h^=h>>16;
		h^=h>>8;
		h^=h>>4;
		h&=7;
		return h;
	}
	float		evaluate(double x, double y)
	{
		static const int log_g2=3,//grid of 8x8 blocks
			g2=1<<log_g2, mask_g2=g2-1,
			g_mask=g_side-1;
		int ax=(int)floor(x), ay=(int)floor(y);
		int hx0=hash3(ax), hx1=hash3(ax+1), hy0=hash3(ay)<<log_g, hy1=hash3(ay+1)<<log_g;
		vec2
			&n_NW=normals[hy1|hx0], &n_NE=normals[hy1|hx1],
			&n_SW=normals[hy0|hx0], &n_SE=normals[hy0|hx1];
		float dx=float(x-ax), dy=float(y-ay);//coordinates inside test cube
		vec2
			v_NW(dx, 1-dy), v_NE(1-dx, 1-dy),
			v_SW(dx, dy), v_SE(1-dx, dy);
		float
			d_NW=v_NW.dot(n_NW), d_NE=v_NE.dot(n_NE),
			d_SW=v_SW.dot(n_SW), d_SE=v_SE.dot(n_SE);
		float
			d_N=mix(d_NW, d_NE, dx),
			d_S=mix(d_SW, d_SE, dx);
		return mix(d_S, d_N, dy);
	}
	float		evaluate(int x, int y)
	{
		static const int log_g2=3,//grid of 8x8 blocks
			g2=1<<log_g2, mask_g2=g2-1,
			g_mask=g_side-1;
		int ax=x>>log_g2, ay=y>>log_g2;//aligned to grid of size 8 blocks
		int hx0=hash3(ax), hx1=hash3(ax+1), hy0=hash3(ay)<<log_g, hy1=hash3(ay+1)<<log_g;
	//	int hx0=hash(ax)&g_mask, hx1=hash(ax+1)&g_mask, hy0=(hash(ay)&g_mask)<<log_g, hy1=(hash(ay+1)&g_mask)<<log_g;
		vec2
			&n_NW=normals[hy1|hx0], &n_NE=normals[hy1|hx1],
			&n_SW=normals[hy0|hx0], &n_SE=normals[hy0|hx1];
		const float inv_side=1.f/g2;
		float dx=(x&mask_g2)*inv_side, dy=(y&mask_g2)*inv_side;//coordinates inside test cube
		vec2
			v_NW(dx, 1-dy), v_NE(1-dx, 1-dy),
			v_SW(dx, dy), v_SE(1-dx, dy);
		float
			d_NW=v_NW.dot(n_NW), d_NE=v_NE.dot(n_NE),
			d_SW=v_SW.dot(n_SW), d_SE=v_SE.dot(n_SE);
		float
			d_N=mix(d_NW, d_NE, dx),
			d_S=mix(d_SW, d_SE, dx);
		return mix(d_S, d_N, dy);
	}
};//pdf_surface;
struct			PerlinNoise3D//range within +-sqrt(3)/2=+-0.867
{
	static const int log_g=3, g_side=1<<log_g, g_size=1<<(log_g+log_g+log_g);
	vec3		normals[g_size];
	std::hash<int> hash;
	void generate_normals()
	{
		const float inv_rmax=2.f/0xFFFFFF;
		for(int k=0;k<g_size;++k)
		{
			auto &n=normals[k];
			do
			{
				n.x=rn24()*inv_rmax-1;
				n.y=rn24()*inv_rmax-1;
				n.z=rn24()*inv_rmax-1;
			}
			while(n.mag_sq()>1);
			n*=1.f/n.magnitude();
		}
	}
	void		seed(unsigned seed)
	{
		srand(seed);
		xoshiro128_seed(rand()<<15|rand(), rand()<<15|rand(), rand()<<15|rand(), rand()<<15|rand());
		generate_normals();
	}
	void		seed(unsigned long long seed)
	{
		srand((unsigned)seed);
		unsigned long long s0=(long long)rand()<<45|(long long)rand()<<30|rand()<<15|rand();
		srand((unsigned)(seed>>32));
		unsigned long long s1=(long long)rand()<<45|(long long)rand()<<30|rand()<<15|rand();
		xoshiro128_seed(s0, s1);
		generate_normals();
	}
	float		evaluate(int x, int y, int z)
	{
		static const int log_g2=3,
			g2=1<<log_g2, mask_g2=g2-1;
		int ax=x>>log_g2, ay=y>>log_g2, az=z>>log_g2;//aligned to grid of size 8 blocks
		static const int g_mask=g_side-1, log_sqg=log_g<<1;
		int hx0=hash(ax		)&g_mask, hy0=(hash(ay	)&g_mask)<<log_g, hz0=(hash(az	)&g_mask)<<log_sqg,
			hx1=hash(ax+1	)&g_mask, hy1=(hash(ay+1)&g_mask)<<log_g, hz1=(hash(az+1)&g_mask)<<log_sqg;
		vec3
			&n_UNW=normals[hz1|hy1|hx0], &n_UNE=normals[hz1|hy1|hx1],
			&n_USW=normals[hz1|hy0|hx0], &n_USE=normals[hz1|hy0|hx1],

			&n_DNW=normals[hz0|hy1|hx0], &n_DNE=normals[hz0|hy1|hx1],
			&n_DSW=normals[hz0|hy0|hx0], &n_DSE=normals[hz0|hy0|hx1];
		const float inv_side=1.f/g2;
		float dx=(x&mask_g2)*inv_side, dy=(y&mask_g2)*inv_side, dz=(z&mask_g2)*inv_side;//coordinates inside test cube
		vec3
			v_UNW(dx, 1-dy, 1-dz), v_UNE(1-dx, 1-dy, 1-dz),
			v_USW(dx, dy, 1-dz), v_USE(1-dx, dy, 1-dz),

			v_DNW(dx, 1-dy, dz), v_DNE(1-dx, 1-dy, dz),
			v_DSW(dx, dy, dz), v_DSE(1-dx, dy, dz);
		float
			d_UNW=v_UNW.dot(n_UNW), d_UNE=v_UNE.dot(n_UNE),
			d_USW=v_USW.dot(n_USW), d_USE=v_USE.dot(n_USE),

			d_DNW=v_DNW.dot(n_DNW), d_DNE=v_DNE.dot(n_DNE),
			d_DSW=v_DSW.dot(n_DSW), d_DSE=v_DSE.dot(n_DSE);
		float
			d_UN=mix(d_UNW, d_UNE, dx),
			d_US=mix(d_USW, d_USE, dx),
			d_DN=mix(d_DNW, d_DNE, dx),
			d_DS=mix(d_DSW, d_DSE, dx),

			d_U=mix(d_US, d_UN, dy),
			d_D=mix(d_DS, d_DN, dy);
		return mix(d_D, d_U, dz);
	}
} pdf_caves;
#undef		mix
#endif
/*struct			BB_PDF
{
	static const int ncomp=4, nterms=ncomp*3;//3 dimenstions
	float amp[nterms], freq[nterms], phase[nterms];
	void seed(unsigned seed)
	{
		srand(seed);
		const float multiplier=0.1f/RAND_MAX, _2pi_rmax=_2pi/RAND_MAX;
		for(int k=0;k<ncomp;++k)
		{
			freq[k]=0.1f+rand()*(k+1)*multiplier;
			amp[k]=1/freq[k];
			phase[k]=rand()*_2pi_rmax;
		}
		for(int k=0;k<ncomp;++k)
		{
			freq[k+4]=0.1f+rand()*(k+1)*multiplier;
			amp[k+4]=1/freq[k+4];
			phase[k+4]=rand()*_2pi_rmax;
		}
		for(int k=0;k<ncomp;++k)
		{
			freq[k+8]=0.1f+rand()*(k+1)*multiplier;
			amp[k+8]=1/freq[k+8];
			phase[k+8]=rand()*_2pi_rmax;
		}
	}
	float evaluate(int x, int y, int z)
	{
		float result=0;
		for(int k=0;k<ncomp;++k)
			result+=amp[k  ]*cos(freq[k  ]*x+phase[k  ]);
		for(int k=0;k<ncomp;++k)
			result+=amp[k+4]*cos(freq[k+4]*y+phase[k+4]);
		for(int k=0;k<ncomp;++k)
			result+=amp[k+8]*cos(freq[k+8]*z+phase[k+8]);
		return result;
	}
} pdf_air, pdf_dirt;//*/
const int		log_cx=4, log_cy=4, log_cz=8,//2,2,8:4x4x256	4,4,8: 256x16x16 blocks per chunk		log_cz should be 8 because of terrain generation
				cx=1<<log_cx, cy=1<<log_cy, cz=1<<log_cz,
				c_size=1<<(log_cx+log_cy+log_cz);
struct			SurfaceInfo
{
	unsigned short idx;
	char DUSNWE;//down MSB, east LSB, true: exposed
};
//void			initialize_slice(char *c, int zstart, int zend, char *pdf, int pdf_size)
//{
//	const int slice_size=cx*cy;
//	for(int kz=zstart;kz<zend;++kz)
//		for(int k=0;k<slice_size;++k)
//			c[slice_size*kz+k]=pdf[rand()%pdf_size];
//}
struct			BBQuadInfo
{
	int count;
	unsigned tx_id;
	BBQuadInfo(){}
	BBQuadInfo(unsigned tx_id):tx_id(tx_id), count(0){}
};
struct			Chunk//65552 bytes
{
	int			x, y;//world coords in blocks
	std::vector<BBQuadInfo> glDrawInfo;
	int			gl_transparent_start_idx;
	unsigned	gl_buffer;
	union
	{
	//	unsigned short c[c_size], d[cz][cx][cy];//chunk blocks: {hichar: light, lochar: type}
		byte c[c_size], d[cz][cx][cy];//chunk blocks
	};
#ifdef LIGHT_SYSTEM
	union
	{
		byte l[c_size], dl[cz][cx][cy];//sky light info: MSB TTTTSSSS LSB, T: torch light, S: skylight
	};
#endif
	std::list<SurfaceInfo> surfaces;
	std::vector<unsigned short> torches;
	void initialize(int x, int y)//function should depend on seed and coordinates
	{
		torches.clear();
		xoshiro128_seed(world_seed, (long long)y<<32|x);
		this->x=x, this->y=y;
		const int log_slice=log_cy+log_cx, slice_size=1<<log_slice;
		memset(c, B_BEDROCK, slice_size);
		for(int kz=1;kz<150;++kz)
	//	for(int kz=1;kz<128;++kz)
		{
			for(int ky=0;ky<cy;++ky)
			{
				for(int kx=0;kx<cx;++kx)
				{
					int idx=kz<<log_slice|ky<<log_cx|kx;
					auto &ck=c[idx];
					float height=float(130+3*pdf_surface.evaluate(0.05*(x+kx), 0.05*(y+ky))+25*pdf_surface.evaluate(0.01*(x+kx), 0.01*(y+ky)));
				//	float height=float(130+5*pdf_surface.evaluate(0.05*(x+kx), 0.05*(y+ky)));//opensimplex noise			1077ms
					if(kz>height+2)
						ck=B_AIR;
					else if(kz>height)
					{
					//	float h2=pdf_surface.evaluate(10*(x+kx), 10*(y+ky));
					//	if(kx>h2)
						int LOL_1=rn24()&0xFF;
						if(LOL_1<250)
							ck=B_AIR;
						else
							ck=B_GLASS;
					}
					else
					{
						float val=pdf_caves.evaluate(x+kx, y+ky, kz);//perlin noise
						if(val>0.2f)
							ck=B_AIR;
						else if(val>0.198f)
							ck=B_TORCH, torches.push_back(idx);
						else if(val>0.19f)
							ck=B_WALL;
						else
							ck=B_DIRT;
						//if(val>0.21f)
						//	ck=B_AIR;
						//else if(val>0.208f)
						//	ck=B_TORCH;
						//else if(val>0.2f)
						//	ck=B_WALL;
						//else
						//	ck=B_DIRT;
					}
				}
			}
		}

		memset(c+slice_size*150, B_AIR, (256-150)*slice_size);
	}
	void update_surfaces(Chunk const *Ec, Chunk const *Wc, Chunk const *Nc, Chunk const *Sc)//called when block(s) change
	{
		SurfaceInfo si;
		surfaces.clear();
		const int log_slice=log_cy+log_cx;
		for(int kz=0;kz<cz;++kz)
		{
			for(int ky=0;ky<cy;++ky)
			{
				for(int kx=0;kx<cx;++kx)
				{
					si.idx=(kz<<log_cy|ky)<<log_cx|kx;
					byte *current=&c[kz<<log_slice|ky<<log_cx|kx];
					if(*current!=B_AIR)
					{
						const byte
							*E=kx+1>=cx	?Ec?&Ec->c[kz<<log_slice|ky		<<log_cx|0]		:nullptr:&c[kz<<log_slice|ky	<<log_cx|(kx+1)],
							*W=kx-1<0	?Wc?&Wc->c[kz<<log_slice|ky		<<log_cx|(cx-1)]:nullptr:&c[kz<<log_slice|ky	<<log_cx|(kx-1)],
							*N=ky+1>=cy	?Nc?&Nc->c[kz<<log_slice|0		<<log_cx|kx]	:nullptr:&c[kz<<log_slice|(ky+1)<<log_cx|kx],
							*S=ky-1<0	?Sc?&Sc->c[kz<<log_slice|(cy-1)	<<log_cx|kx]	:nullptr:&c[kz<<log_slice|(ky-1)<<log_cx|kx],
							*U=kz+1>=cz	?nullptr:&c[(kz+1)<<log_slice|ky<<log_cx|kx],
							*D=kz-1<0	?nullptr:&c[(kz-1)<<log_slice|ky<<log_cx|kx];
						//if(D&&*D==B_GLASS||U&&*U==B_GLASS||S&&*S==B_GLASS||N&&*N==B_GLASS||W&&*W==B_GLASS||E&&*E==B_GLASS)
						//	int LOL_1=0;
						si.DUSNWE=(D&&*D<=B_GLASS)<<5|(!U||*U<=B_GLASS)<<4|(S&&*S<=B_GLASS)<<3|(N&&*N<=B_GLASS)<<2|(W&&*W<=B_GLASS)<<1|(E&&*E<=B_GLASS);
						if(si.DUSNWE)
							surfaces.push_back(si);
					}
				}
			}
		}
	}
};
const int		log_dd=
#ifdef _DEBUG
				4
#else
				7
#endif
				, draw_distance=1<<log_dd,//in blocks		2: 4|	4: 16 (min)		5: 32	6: 64	7: 128		BROKEN: 8: 256, 10: 1024	//TODO: adjustable draw distance
				log_level_dx=log_dd+1-log_cx, log_level_dy=log_dd+1-log_cy,
				level_dx=draw_distance<<1>>log_cx, level_dy=draw_distance<<1>>log_cy, level_size=level_dx*level_dy;//in chunks
std::vector<Chunk> chunks;
int				VX=0, VY=0;
//int				arrayXstart=-draw_distance, arrayYstart=-draw_distance;
inline void		update_surfaces(int gx, int gy)
{
	chunks[gy<<log_level_dx|gx].update_surfaces(
		gx+1<level_dx?&chunks[gy<<log_level_dx|(gx+1)]:nullptr, gx-1>=0?&chunks[gy<<log_level_dx|(gx-1)]:nullptr,
		gy+1<level_dy?&chunks[(gy+1)<<log_level_dx|gx]:nullptr, gy-1>=0?&chunks[(gy-1)<<log_level_dx|gx]:nullptr);
}
inline byte*	get(int wx, int wy, int z)
{
	int ax=wx+draw_distance, ay=wy+draw_distance;//array coordinates
	int gx=ax>>log_cx, gy=ay>>log_cy;
	if(gx>=0&&gx<level_dx&&gy>=0&&gy<level_dy&&z>=0&&z<cz)
	{
		const int cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cy+log_cx;
		int lx=ax&cx_mask, ly=ay&cy_mask;
		auto &ch=chunks[gy<<log_level_dx|gx];
		return ch.c+(z<<log_slice|ly<<log_cx|lx);
	}
	return nullptr;
}
inline void		get(int ax, int ay, int z, byte *&block, byte *&lightinfo)//array coordinates
{
	int gx=ax>>log_cx, gy=ay>>log_cy;
	if(gx>=0&&gx<level_dx&&gy>=0&&gy<level_dy&&z>=0&&z<cz)
	{
		const int cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cy+log_cx;
		int lx=ax&cx_mask, ly=ay&cy_mask;
		auto &ch=chunks[gy<<log_level_dx|gx];
		int idx=z<<log_slice|ly<<log_cx|lx;
		block=ch.c+idx, lightinfo=ch.l+idx;
	}
}
inline void		get_world(int wx, int wy, int z, byte *&block, byte *&lightinfo)
{
	int ax=wx+draw_distance, ay=wy+draw_distance;//array coordinates
	int gx=ax>>log_cx, gy=ay>>log_cy;
	if(gx>=0&&gx<level_dx&&gy>=0&&gy<level_dy&&z>=0&&z<cz)
	{
		const int cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cy+log_cx;
		int lx=ax&cx_mask, ly=ay&cy_mask;
		auto &ch=chunks[gy<<log_level_dx|gx];
		int idx=z<<log_slice|ly<<log_cx|lx;
		block=ch.c+idx, lightinfo=ch.l+idx;
	}
}
inline int		get(int wx, int wy, int z, int &ch_x, int &ch_y, int &b_idx, int &lx, int &ly)//world coordinates, returns chunk idx -1 if OOB
{
	int ax=wx+draw_distance, ay=wy+draw_distance;//array coordinates
	int gx=ax>>log_cx, gy=ay>>log_cy;
	if(gx>=0&&gx<level_dx&&gy>=0&&gy<level_dy&&z>=0&&z<cz)
	{
		const int cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cy+log_cx;
		lx=ax&cx_mask, ly=ay&cy_mask;
		ch_x=gx, ch_y=gy;
		int ch_idx=gy<<log_level_dx|gx;
		b_idx=z<<log_slice|ly<<log_cx|lx;
		return ch_idx;
	}
	return -1;
}
#if 0
//inline char*	get(int gx, int gy, int lx, int ly, int z)//inter-chunk coordinates, in-chunk coordinates
//{
//	if(gx>=0&&gx<level_dx&&gy>=0&&gy<level_dy&&z>=0&&z<cz)
//	{
//		const int cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cy+log_cx;
//		auto &ch=chunks[gy<<log_level_dx|gx];
//		return ch.c+(z<<log_slice|ly<<log_cx|lx);
//	}
//	return nullptr;
//}
char			light_iterations=0;
void			calculate_lights_rain()
{
	const int
		log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,//total number of blocks in x & y directions
		log_slice=log_cy+log_cx, slice_size=1<<log_slice;
	for(int gy=0;gy<level_dy;++gy)//initiate light in transparent sky blocks
	{
		for(int gx=0;gx<level_dx;++gx)
		{
			auto &ch=chunks[gy<<log_level_dx|gx];
			memset(ch.l, 0, c_size);
			for(int ky=0;ky<cy;++ky)
				for(int kx=0;kx<cx;++kx)//set all sky facing transparent blocks to maximum light
					for(int kz=cz-1;kz>=0&&ch.c[kz<<log_slice|ky<<log_cx|kx]<=B_GLASS;--kz)
						ch.l[kz<<log_slice|ky<<log_cx|kx]=15;
		}
	}
}
void			calculate_lights_iterate()
{
	const int
		log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,//total number of blocks in x & y directions
		log_slice=log_cy+log_cx, slice_size=1<<log_slice;
	for(int gy=0;gy<level_dy;++gy)
	{
		for(int gx=0;gx<level_dx;++gx)
		{
			auto &ch=chunks[gy<<log_level_dx|gx];
			for(int z=0;z<cz;++z)
			{
				for(int ly=0;ly<cy;++ly)
				{
					for(int lx=0;lx<cx;++lx)
					{
						int idx=z<<log_slice|ly<<log_cx|lx;
						if(ch.c[idx]<=B_GLASS)
						{
							byte &Cl=ch.l[idx];//current lightinfo
							int ax=gx<<log_cx|lx, ay=gy<<log_cy|ly;
							byte
								*Ub=nullptr, *Db=nullptr, *Nb=nullptr, *Sb=nullptr, *Eb=nullptr, *Wb=nullptr,
								*Ul=nullptr, *Dl=nullptr, *Nl=nullptr, *Sl=nullptr, *El=nullptr, *Wl=nullptr;
							get(ax, ay, z+1, Ub, Ul), get(ax, ay, z-1, Db, Dl);
							get(ax, ay+1, z, Nb, Nl), get(ax, ay-1, z, Sb, Sl);
							get(ax+1, ay, z, Eb, El), get(ax-1, ay, z, Wb, Wl);
							//if(gx==1&&gy==2&&lx==9&&ly==0&&z==122)//
							//	int LOL_1=0;
							//if(gx==1&&gy==2&&lx==9&&ly==0&&z==121)//
							//	int LOL_1=0;
							//if(gx==1&&gy==2&&lx==7&&ly==4&&z==122)//
							//	int LOL_1=0;
							//auto Cl0=Cl;//
							if(Ub&&*Ub<=B_GLASS&&Cl<*Ul-1)
								Cl=*Ul-1;
							if(Db&&*Db<=B_GLASS&&Cl<*Dl-1)
								Cl=*Dl-1;
							if(Nb&&*Nb<=B_GLASS&&Cl<*Nl-1)
								Cl=*Nl-1;
							if(Sb&&*Sb<=B_GLASS&&Cl<*Sl-1)
								Cl=*Sl-1;
							if(Eb&&*Eb<=B_GLASS&&Cl<*El-1)
								Cl=*El-1;
							if(Wb&&*Wb<=B_GLASS&&Cl<*Wl-1)
								Cl=*Wl-1;
							//if(Cl>Cl0)//
							//	int LOL_1=0;//
						//	Cl=(lx+light_iterations)&15;//
						}
					}
				}
			}
		}
	}
}
#endif
int				calculate_lights_timeout=25000000;//25M: 26fps, 55M: 20fps, 117M: 15fps
typedef std::pair<unsigned short, unsigned short> CB;//chunk, block
std::queue<CB>	q;
bool			lights_done=false;
int				start_gy=0, start_gx=0, start_ly=0, start_lx=0;
void			update_lights(unsigned short ch_idx, unsigned short b_idx)
{
	const int
		log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,//total number of blocks in x & y directions
		log_slice=log_cy+log_cx, slice_size=1<<log_slice,
		cx_mask=cx-1, cy_mask=cy-1,
		level_dx_mask=level_dx-1, level_dy_mask=level_dy-1;
	std::queue<CB> q2;//recalculate lights in the neighborhood
	q2.push(CB(ch_idx, b_idx));
	while(q2.size())//breadth-first traversal - torch lights
	{
		auto &front=q2.front();
		unsigned short ch2_idx=front.first, gx2=ch2_idx&level_dx_mask, gy2=ch2_idx>>log_level_dx&level_dy_mask;
		Chunk *ch2=&chunks[ch2_idx];
		int current_idx=front.second, x2=current_idx&cx_mask, y2=current_idx>>log_cx&cy_mask, z2=current_idx>>log_slice;//in-chunk coordinates
		q2.pop();
		int flood=ch2->l[z2<<log_slice|y2<<log_cx|x2]>>4, flood_m1=flood-1;
		unsigned short idx;
		if(x2-1>=0)							//WEST: -x
		{
			idx=z2<<log_slice|y2<<log_cx|(x2-1);
			if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
				ch2->l[idx]=flood_m1<<4|ch2->l[idx]&15, q2.push(CB(ch2_idx, idx));
		}
		else if(gx2-1>=0)
		{
			unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
			auto ch3=&chunks[ch3_idx];
			idx=z2<<log_slice|y2<<log_cx|cx_mask;
			if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
				ch3->l[idx]=flood_m1<<4|ch3->l[idx]&15, q2.push(CB(ch3_idx, idx));
		}
		if(x2+1<cx)							//EAST: +x
		{
			idx=z2<<log_slice|y2<<log_cx|(x2+1);
			if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
				ch2->l[idx]=flood_m1<<4|ch2->l[idx]&15, q2.push(CB(ch2_idx, idx));
		}
		else if(gx2+1<level_dx)
		{
			unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
			auto ch3=&chunks[ch3_idx];
			idx=z2<<log_slice|y2<<log_cx|0;
			if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
				ch3->l[idx]=flood_m1<<4|ch3->l[idx]&15, q2.push(CB(ch3_idx, idx));
		}
		if(y2-1>=0)							//SOUTH: -y
		{
			idx=z2<<log_slice|(y2-1)<<log_cx|x2;
			if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
				ch2->l[idx]=flood_m1<<4|ch2->l[idx]&15, q2.push(CB(ch2_idx, idx));
		}
		else if(gy2-1>=0)
		{
			unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
			auto ch3=&chunks[ch3_idx];
			idx=z2<<log_slice|cy_mask<<log_cx|x2;
			if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
				ch3->l[idx]=flood_m1<<4|ch3->l[idx]&15, q2.push(CB(ch3_idx, idx));
		}
		if(y2+1<cy)							//NORTH: +y
		{
			idx=z2<<log_slice|(y2+1)<<log_cx|x2;
			if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
				ch2->l[idx]=flood_m1<<4|ch2->l[idx]&15, q2.push(CB(ch2_idx, idx));
		}
		else if(gy2+1<level_dy)
		{
			unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
			auto ch3=&chunks[ch3_idx];
			idx=z2<<log_slice|0<<log_cx|x2;
			if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
				ch3->l[idx]=flood_m1<<4|ch3->l[idx]&15, q2.push(CB(ch3_idx, idx));
		}
		if(z2-1>=0)							//DOWN: -z
		{
			idx=(z2-1)<<log_slice|y2<<log_cx|x2;
			if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
				ch2->l[idx]=flood_m1<<4|ch2->l[idx]&15, q2.push(CB(ch2_idx, idx));
		}
		if(z2+1<cz)							//UP: +z
		{
			idx=(z2+1)<<log_slice|y2<<log_cx|x2;
			if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
				ch2->l[idx]=flood_m1<<4|ch2->l[idx]&15, q2.push(CB(ch2_idx, idx));
		}
	}
	q2.push(CB(ch_idx, b_idx));
	while(q2.size())//breadth-first traversal - skylight
	{
		auto &front=q2.front();
		unsigned short ch2_idx=front.first, gx2=ch2_idx&level_dx_mask, gy2=ch2_idx>>log_level_dx&level_dy_mask;
		Chunk *ch2=&chunks[ch2_idx];
		int current_idx=front.second, x2=current_idx&cx_mask, y2=current_idx>>log_cx&cy_mask, z2=current_idx>>log_slice;//in-chunk coordinates
		q2.pop();
		int flood=ch2->l[z2<<log_slice|y2<<log_cx|x2]&15, flood_m1=flood-1;
		unsigned short idx;
		if(flood==15)
		{
			for(int z3=z2;z3>=0;--z3)
			{
				idx=z3<<log_slice|y2<<log_cx|x2;
				if(ch2->c[idx]>B_GLASS)
					break;
				ch2->l[idx]|=15;
				if(x2-1>=0)							//WEST: -x
				{
					idx=z3<<log_slice|y2<<log_cx|(x2-1);
					if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
						ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
				}
				else if(gx2-1>=0)
				{
					unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
					auto ch3=&chunks[ch3_idx];
					idx=z3<<log_slice|y2<<log_cx|cx_mask;
					if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
						ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
				}
				if(x2+1<cx)							//EAST: +x
				{
					idx=z3<<log_slice|y2<<log_cx|(x2+1);
					if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
						ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
				}
				else if(gx2+1<level_dx)
				{
					unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
					auto ch3=&chunks[ch3_idx];
					idx=z3<<log_slice|y2<<log_cx|0;
					if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
						ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
				}
				if(y2-1>=0)							//SOUTH: -y
				{
					idx=z3<<log_slice|(y2-1)<<log_cx|x2;
					if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
						ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
				}
				else if(gy2-1>=0)
				{
					unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
					auto ch3=&chunks[ch3_idx];
					idx=z3<<log_slice|cy_mask<<log_cx|x2;
					if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
						ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
				}
				if(y2+1<cy)							//NORTH: +y
				{
					idx=z3<<log_slice|(y2+1)<<log_cx|x2;
					if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
						ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
				}
				else if(gy2+1<level_dy)
				{
					unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
					auto ch3=&chunks[ch3_idx];
					idx=z3<<log_slice|0<<log_cx|x2;
					if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
						ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
				}
			}
		}
		else//flood<15
		{
			if(x2-1>=0)							//WEST: -x
			{
				idx=z2<<log_slice|y2<<log_cx|(x2-1);
				if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
					ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
			}
			else if(gx2-1>=0)
			{
				unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
				auto ch3=&chunks[ch3_idx];
				idx=z2<<log_slice|y2<<log_cx|cx_mask;
				if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
					ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
			}
			if(x2+1<cx)							//EAST: +x
			{
				idx=z2<<log_slice|y2<<log_cx|(x2+1);
				if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
					ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
			}
			else if(gx2+1<level_dx)
			{
				unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
				auto ch3=&chunks[ch3_idx];
				idx=z2<<log_slice|y2<<log_cx|0;
				if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
					ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
			}
			if(y2-1>=0)							//SOUTH: -y
			{
				idx=z2<<log_slice|(y2-1)<<log_cx|x2;
				if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
					ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
			}
			else if(gy2-1>=0)
			{
				unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
				auto ch3=&chunks[ch3_idx];
				idx=z2<<log_slice|cy_mask<<log_cx|x2;
				if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
					ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
			}
			if(y2+1<cy)							//NORTH: +y
			{
				idx=z2<<log_slice|(y2+1)<<log_cx|x2;
				if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
					ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
			}
			else if(gy2+1<level_dy)
			{
				unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
				auto ch3=&chunks[ch3_idx];
				idx=z2<<log_slice|0<<log_cx|x2;
				if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
					ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q2.push(CB(ch3_idx, idx));
			}
			if(z2-1>=0)							//DOWN: -z
			{
				idx=(z2-1)<<log_slice|y2<<log_cx|x2;
				if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
					ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
			}
			if(z2+1<cz)							//UP: +z
			{
				idx=(z2+1)<<log_slice|y2<<log_cx|x2;
				if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
					ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q2.push(CB(ch2_idx, idx));
			}
		}
	}
}
void			calculate_lights_v3(bool start_over)
{
	long long t_start=__rdtsc();
	const int
		log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,//total number of blocks in x & y directions
		log_slice=log_cy+log_cx, slice_size=1<<log_slice,
		cx_mask=cx-1, cy_mask=cy-1,
		level_dx_mask=level_dx-1, level_dy_mask=level_dy-1;
	bool fresh_entry=!start_over;
	if(start_over)
	{
		for(int k=0, kEnd=chunks.size();k<kEnd;++k)//reset all lightinfo to zero
			memset(chunks[k].l, 0, c_size);
		//for(int gy=0;gy<level_dy;++gy)
		//	for(int gx=0;gx<level_dx;++gx)
		//		memset(chunks[gy<<log_level_dx|gx].l, 0, c_size);
		std::queue<CB> q2;
		for(int k=0, kEnd=chunks.size();k<kEnd;++k)//calculate all torch lights
		{
			auto &ch=chunks[k];
			for(int kt=0, ktEnd=ch.torches.size();kt<ktEnd;++kt)
			{
				auto t_idx=ch.torches[kt];
				ch.l[t_idx]=15<<4;
				q2.push(CB(k, t_idx));
				while(q2.size())//breadth-first traversal
				{
					auto &front=q2.front();
					unsigned short ch2_idx=front.first, gx2=ch2_idx&level_dx_mask, gy2=ch2_idx>>log_level_dx&level_dy_mask;
					Chunk *ch2=&chunks[ch2_idx];
					int current_idx=front.second, x2=current_idx&cx_mask, y2=current_idx>>log_cx&cy_mask, z2=current_idx>>log_slice;//in-chunk coordinates
					q2.pop();
					int flood=ch2->l[z2<<log_slice|y2<<log_cx|x2]>>4, flood_m1=flood-1;
					unsigned short idx;
					if(x2-1>=0)							//WEST: -x
					{
						idx=z2<<log_slice|y2<<log_cx|(x2-1);
						if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
							ch2->l[idx]=flood_m1<<4, q2.push(CB(ch2_idx, idx));
					}
					else if(gx2-1>=0)
					{
						unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
						auto ch3=&chunks[ch3_idx];
						idx=z2<<log_slice|y2<<log_cx|cx_mask;
						if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
							ch3->l[idx]=flood_m1<<4, q2.push(CB(ch3_idx, idx));
					}
					if(x2+1<cx)							//EAST: +x
					{
						idx=z2<<log_slice|y2<<log_cx|(x2+1);
						if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
							ch2->l[idx]=flood_m1<<4, q2.push(CB(ch2_idx, idx));
					}
					else if(gx2+1<level_dx)
					{
						unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
						auto ch3=&chunks[ch3_idx];
						idx=z2<<log_slice|y2<<log_cx|0;
						if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
							ch3->l[idx]=flood_m1<<4, q2.push(CB(ch3_idx, idx));
					}
					if(y2-1>=0)							//SOUTH: -y
					{
						idx=z2<<log_slice|(y2-1)<<log_cx|x2;
						if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
							ch2->l[idx]=flood_m1<<4, q2.push(CB(ch2_idx, idx));
					}
					else if(gy2-1>=0)
					{
						unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
						auto ch3=&chunks[ch3_idx];
						idx=z2<<log_slice|cy_mask<<log_cx|x2;
						if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
							ch3->l[idx]=flood_m1<<4, q2.push(CB(ch3_idx, idx));
					}
					if(y2+1<cy)							//NORTH: +y
					{
						idx=z2<<log_slice|(y2+1)<<log_cx|x2;
						if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
							ch2->l[idx]=flood_m1<<4, q2.push(CB(ch2_idx, idx));
					}
					else if(gy2+1<level_dy)
					{
						unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
						auto ch3=&chunks[ch3_idx];
						idx=z2<<log_slice|0<<log_cx|x2;
						if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]>>4)<flood_m1)
							ch3->l[idx]=flood_m1<<4, q2.push(CB(ch3_idx, idx));
					}
					if(z2-1>=0)							//DOWN: -z
					{
						idx=(z2-1)<<log_slice|y2<<log_cx|x2;
						if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
							ch2->l[idx]=flood_m1<<4, q2.push(CB(ch2_idx, idx));
					}
					if(z2+1<cz)							//UP: +z
					{
						idx=(z2+1)<<log_slice|y2<<log_cx|x2;
						if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]>>4)<flood_m1)
							ch2->l[idx]=flood_m1<<4, q2.push(CB(ch2_idx, idx));
					}
				}
			}
		}
		q=std::queue<CB>();
		start_gy=0, start_gx=0, start_ly=0, start_lx=0;
	}
	for(int gy=start_gy;gy<level_dy;++gy)//calculate light on block faces
	{
		start_gy=0;
		for(int gx=start_gx;gx<level_dx;++gx)
		{
			start_gx=0;
			unsigned short ch_idx=gy<<log_level_dx|gx;
			auto &ch=chunks[ch_idx];
			for(int ky=start_ly;ky<cy;++ky)
			{
				start_ly=0;
				for(int kx=start_lx;kx<cx;++kx)
				{
					start_lx=0;
					int start_idx=(cz-1)<<log_slice|ky<<log_cx|kx;
					if(ch.c[start_idx]<=B_GLASS)//top block is transparent
					{
						ch.l[start_idx]|=15;
						if(!fresh_entry)
							q.push(CB(ch_idx, start_idx));
						fresh_entry=false;
						while(q.size())//breadth-first traversal
						{
							if(__rdtsc()-t_start>calculate_lights_timeout)
							{
								start_gy=gy, start_gx=gx, start_ly=ky, start_lx=kx;
								return;
							}
							auto &front=q.front();
							unsigned short ch2_idx=front.first, gx2=ch2_idx&level_dx_mask, gy2=ch2_idx>>log_level_dx&level_dy_mask;
							Chunk *ch2=&chunks[ch2_idx];
							int current_idx=front.second, x2=current_idx&cx_mask, y2=current_idx>>log_cx&cy_mask, z2=current_idx>>log_slice;//in-chunk coordinates
							q.pop();
							int flood=ch2->l[z2<<log_slice|y2<<log_cx|x2]&15, flood_m1=flood-1;
							unsigned short idx;
							if(flood==15)
							{
								for(int z3=z2;z3>=0;--z3)
								{
									idx=z3<<log_slice|y2<<log_cx|x2;
									if(ch2->c[idx]>B_GLASS)
										break;
									ch2->l[idx]|=15;
									if(x2-1>=0)							//WEST: -x
									{
										idx=z3<<log_slice|y2<<log_cx|(x2-1);
										if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
											ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gx2-1>=0)
									{
										unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|y2<<log_cx|cx_mask;
										if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
											ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
									}
									if(x2+1<cx)							//EAST: +x
									{
										idx=z3<<log_slice|y2<<log_cx|(x2+1);
										if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
											ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gx2+1<level_dx)
									{
										unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|y2<<log_cx|0;
										if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
											ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
									}
									if(y2-1>=0)							//SOUTH: -y
									{
										idx=z3<<log_slice|(y2-1)<<log_cx|x2;
										if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
											ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gy2-1>=0)
									{
										unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|cy_mask<<log_cx|x2;
										if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
											ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
									}
									if(y2+1<cy)							//NORTH: +y
									{
										idx=z3<<log_slice|(y2+1)<<log_cx|x2;
										if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
											ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gy2+1<level_dy)
									{
										unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|0<<log_cx|x2;
										if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
											ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
									}
								}
							}
							else//flood<15
							{
								if(x2-1>=0)							//WEST: -x
								{
									idx=z2<<log_slice|y2<<log_cx|(x2-1);
									if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
										ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gx2-1>=0)
								{
									unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|y2<<log_cx|cx_mask;
									if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
										ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(x2+1<cx)							//EAST: +x
								{
									idx=z2<<log_slice|y2<<log_cx|(x2+1);
									if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
										ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gx2+1<level_dx)
								{
									unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|y2<<log_cx|0;
									if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
										ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(y2-1>=0)							//SOUTH: -y
								{
									idx=z2<<log_slice|(y2-1)<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
										ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gy2-1>=0)
								{
									unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|cy_mask<<log_cx|x2;
									if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
										ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(y2+1<cy)							//NORTH: +y
								{
									idx=z2<<log_slice|(y2+1)<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
										ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gy2+1<level_dy)
								{
									unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|0<<log_cx|x2;
									if(ch3->c[idx]<=B_GLASS&&(ch3->l[idx]&15)<flood_m1)
										ch3->l[idx]=ch3->l[idx]&0xF0|flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(z2-1>=0)							//DOWN: -z
								{
									idx=(z2-1)<<log_slice|y2<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
										ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
								}
								if(z2+1<cz)							//UP: +z
								{
									idx=(z2+1)<<log_slice|y2<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&(ch2->l[idx]&15)<flood_m1)
										ch2->l[idx]=ch2->l[idx]&0xF0|flood_m1, q.push(CB(ch2_idx, idx));
								}
							}
						}
					}
				}
			}
		}
	}
	lights_done=true;
}
#if 0
void			calculate_lights_v2()
{
	const int
		log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,//total number of blocks in x & y directions
		log_slice=log_cy+log_cx, slice_size=1<<log_slice;
	for(int gy=0;gy<level_dy;++gy)//reset all lightinfo to zero
		for(int gx=0;gx<level_dx;++gx)
			memset(chunks[gy<<log_level_dx|gx].l, 0, c_size);
	const int
		cx_mask=cx-1, cy_mask=cy-1,
		level_dx_mask=level_dx-1, level_dy_mask=level_dy-1;
	typedef std::pair<unsigned short, unsigned short> CB;//chunk, block
	std::queue<CB> q;
	for(int gy=0;gy<level_dy;++gy)//calculate light on block faces
	{
		for(int gx=0;gx<level_dx;++gx)
		{
			unsigned short ch_idx=gy<<log_level_dx|gx;
			auto &ch=chunks[ch_idx];
			for(int ky=0;ky<cy;++ky)
			{
				for(int kx=0;kx<cx;++kx)
				{
					if(ch.c[(cz-1)<<log_slice|ky<<log_cx|kx]<=B_GLASS)//top block is transparent
					{
						ch.l[(cz-1)<<log_slice|ky<<log_cx|kx]=15;
						q.push(CB(ch_idx, (cz-1)<<log_slice|ky<<log_cx|kx));
						while(q.size())//breadth-first traversal
						{
							auto &front=q.front();
							unsigned short ch2_idx=front.first, gx2=ch2_idx&level_dx_mask, gy2=ch2_idx>>log_level_dx&level_dy_mask;
							Chunk *ch2=&chunks[ch2_idx];
							int current_idx=front.second, x2=current_idx&cx_mask, y2=current_idx>>log_cx&cy_mask, z2=current_idx>>log_slice;//in-chunk coordinates
							q.pop();
							int flood=ch2->l[z2<<log_slice|y2<<log_cx|x2], flood_m1=flood-1;
							unsigned short idx;
							if(flood==15)
							{
								for(int z3=z2;z3>=0;--z3)
								{
									idx=z3<<log_slice|y2<<log_cx|x2;
									if(ch2->c[idx]>B_GLASS)
										break;
									ch2->l[idx]=15;
									if(x2-1>=0)							//WEST: -x
									{
										idx=z3<<log_slice|y2<<log_cx|(x2-1);
										if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
											ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gx2-1>=0)
									{
										unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|y2<<log_cx|cx_mask;
										if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
											ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
									}
									if(x2+1<cx)							//EAST: +x
									{
										idx=z3<<log_slice|y2<<log_cx|(x2+1);
										if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
											ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gx2+1<level_dx)
									{
										unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|y2<<log_cx|0;
										if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
											ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
									}
									if(y2-1>=0)							//SOUTH: -y
									{
										idx=z3<<log_slice|(y2-1)<<log_cx|x2;
										if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
											ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gy2-1>=0)
									{
										unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|cy_mask<<log_cx|x2;
										if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
											ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
									}
									if(y2+1<cy)							//NORTH: +y
									{
										idx=z3<<log_slice|(y2+1)<<log_cx|x2;
										if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
											ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
									}
									else if(gy2+1<level_dy)
									{
										unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
										auto ch3=&chunks[ch3_idx];
										idx=z3<<log_slice|0<<log_cx|x2;
										if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
											ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
									}
								}
							}
							else//flood<15
							{
								if(x2-1>=0)							//WEST: -x
								{
									idx=z2<<log_slice|y2<<log_cx|(x2-1);
									if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
										ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gx2-1>=0)
								{
									unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|y2<<log_cx|cx_mask;
									if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
										ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(x2+1<cx)							//EAST: +x
								{
									idx=z2<<log_slice|y2<<log_cx|(x2+1);
									if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
										ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gx2+1<level_dx)
								{
									unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|y2<<log_cx|0;
									if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
										ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(y2-1>=0)							//SOUTH: -y
								{
									idx=z2<<log_slice|(y2-1)<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
										ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gy2-1>=0)
								{
									unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|cy_mask<<log_cx|x2;
									if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
										ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(y2+1<cy)							//NORTH: +y
								{
									idx=z2<<log_slice|(y2+1)<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
										ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
								}
								else if(gy2+1<level_dy)
								{
									unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
									auto ch3=&chunks[ch3_idx];
									idx=z2<<log_slice|0<<log_cx|x2;
									if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
										ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
								}
								if(z2-1>=0)							//DOWN: -z
								{
									idx=(z2-1)<<log_slice|y2<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
										ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
								}
								if(z2+1<cz)							//UP: +z
								{
									idx=(z2+1)<<log_slice|y2<<log_cx|x2;
									if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
										ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
								}
							}
						}
					}
				}
			}
		}
	}
}
void			calculate_lights()
{
	const int
		log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,//total number of blocks in x & y directions
		log_slice=log_cy+log_cx, slice_size=1<<log_slice;
	//unsigned char z_ground[1<<(log_by+log_bx)];
	for(int gy=0;gy<level_dy;++gy)//reset all lightinfo to zero
		for(int gx=0;gx<level_dx;++gx)
			memset(chunks[gy<<log_level_dx|gx].l, 0, c_size);
	//for(int gy=0;gy<level_dy;++gy)//initiate light in transparent sky blocks
	//{
	//	for(int gx=0;gx<level_dx;++gx)
	//	{
	//		auto &ch=chunks[gy<<log_level_dx|gx];
	//		memset(ch.l, 0, c_size);
	//		for(int ky=0;ky<cy;++ky)
	//		{
	//			for(int kx=0;kx<cx;++kx)//find ground elevation
	//			{
	//				int kz=cz-1;
	//				for(;kz>=0&&ch.d[kz][ky][kx]<=B_GLASS;--kz);//skip all transparent blocks
	//				z_ground[(gy<<log_cy|ky)<<log_bx|(gx<<log_cx|kx)]=kz;
	//			}
	//		}
	//	}
	//}
	const int
		cx_mask=cx-1, cy_mask=cy-1,
		level_dx_mask=level_dx-1, level_dy_mask=level_dy-1;
	typedef std::pair<unsigned short, unsigned short> CB;//chunk, block
	std::queue<CB> q;
	for(int gy=0;gy<level_dy;++gy)//calculate light on block faces
	{
		for(int gx=0;gx<level_dx;++gx)
		{
			unsigned short ch_idx=gy<<log_level_dx|gx;
			auto &ch=chunks[ch_idx];
			for(int ky=0;ky<cy;++ky)
			{
				for(int kx=0;kx<cx;++kx)
				{
				//	auto current_z_a=z_ground[(gy<<log_cy|ky)<<log_bx|(gx<<log_cx|kx)]+1;
					if(ch.c[(cz-1)<<log_slice|ky<<log_cx|kx]<=B_GLASS)
				//	if(current_z_a<cz)//transparent block just above opaque ground is below cz
					{
						ch.l[(cz-1)<<log_slice|ky<<log_cx|kx]=15;
						q.push(CB(ch_idx, (cz-1)<<log_slice|ky<<log_cx|kx));
					//	q.push(CB(ch_idx, current_z_g<<log_slice|ky<<log_cx|kx));
						while(q.size())//breadth-first traversal
						{
							auto &front=q.front();
							unsigned short ch2_idx=front.first, gx2=ch2_idx&level_dx_mask, gy2=ch2_idx>>log_level_dx&level_dy_mask;
							Chunk *ch2=&chunks[ch2_idx];
							int current_idx=front.second, x2=current_idx&cx_mask, y2=current_idx>>log_cx&cy_mask, z2=current_idx>>log_slice;//in-chunk coordinates
							q.pop();
							int flood=ch2->l[z2<<log_slice|y2<<log_cx|x2], flood_m1=flood-1, flood_down=flood-(flood<15);
							unsigned short idx;
							if(x2-1>=0)							//WEST: -x
							{
								idx=z2<<log_slice|y2<<log_cx|(x2-1);
								if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
									ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
							}
							else if(gx2-1>=0)
							{
								unsigned short ch3_idx=gy2<<log_level_dx|(gx2-1);
								auto ch3=&chunks[ch3_idx];
								idx=z2<<log_slice|y2<<log_cx|cx_mask;
								if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
									ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
							}
							if(x2+1<cx)							//EAST: +x
							{
								idx=z2<<log_slice|y2<<log_cx|(x2+1);
								if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
									ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
							}
							else if(gx2+1<level_dx)
							{
								unsigned short ch3_idx=gy2<<log_level_dx|(gx2+1);
								auto ch3=&chunks[ch3_idx];
								idx=z2<<log_slice|y2<<log_cx|0;
								if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
									ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
							}
							if(y2-1>=0)							//SOUTH: -y
							{
								idx=z2<<log_slice|(y2-1)<<log_cx|x2;
								if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
									ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
							}
							else if(gy2-1>=0)
							{
								unsigned short ch3_idx=(gy2-1)<<log_level_dx|gx2;
								auto ch3=&chunks[ch3_idx];
								idx=z2<<log_slice|cy_mask<<log_cx|x2;
								if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
									ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
							}
							if(y2+1<cy)							//NORTH: +y
							{
								idx=z2<<log_slice|(y2+1)<<log_cx|x2;
								if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
									ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
							}
							else if(gy2+1<level_dy)
							{
								unsigned short ch3_idx=(gy2+1)<<log_level_dx|gx2;
								auto ch3=&chunks[ch3_idx];
								idx=z2<<log_slice|0<<log_cx|x2;
								if(ch3->c[idx]<=B_GLASS&&ch3->l[idx]<flood_m1)
									ch3->l[idx]=flood_m1, q.push(CB(ch3_idx, idx));
							}
							if(z2-1>=0)							//DOWN: -z
							{
								idx=(z2-1)<<log_slice|y2<<log_cx|x2;
								if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_down)
									ch2->l[idx]=flood_down, q.push(CB(ch2_idx, idx));
							}
							if(z2+1<cz)							//UP: +z
							{
								idx=(z2+1)<<log_slice|y2<<log_cx|x2;
								if(ch2->c[idx]<=B_GLASS&&ch2->l[idx]<flood_m1)
									ch2->l[idx]=flood_m1, q.push(CB(ch2_idx, idx));
							}
						}
					}
				}
			}
		}
	}
}
#endif
void			explore(int x, int y, int dx, int dy)//coords & displacement in blocks
{
	if(chunks.size()!=level_size)//initiation
	{
		for(int k=0, kEnd=chunks.size();k<kEnd;++k)//
			glDeleteBuffers(1, &chunks[k].gl_buffer);//
		chunks.resize(level_size);
		for(int k=0, kEnd=chunks.size();k<kEnd;++k)
			glGenBuffers(1, &chunks[k].gl_buffer);

		int xstart=VX-draw_distance, ystart=VY-draw_distance;
		xstart-=xstart&cx-1;
		ystart-=ystart&cy-1;
		//arrayXstart=-draw_distance, arrayYstart=-draw_distance;//start at center
		//arrayXstart-=arrayXstart&cx_mask;
		//arrayYstart-=arrayYstart&cy_mask;
		//int xstart=x-draw_distance, ystart=y-draw_distance;
		//xstart-=xstart&cx-1;
		//ystart-=ystart&cy-1;
		prof_add("Gen. ch. start");
		for(int ky=0;ky<level_dy;++ky)//generate all active chunks
			for(int kx=0;kx<level_dx;++kx)
				chunks[ky<<log_level_dx|kx].initialize(xstart+(kx<<log_cx), ystart+(ky<<log_cy));
		prof_add("Generate chunks");
		for(int ky=0;ky<level_dy;++ky)//cache all exposed quads
			for(int kx=0;kx<level_dx;++kx)
				update_surfaces(kx, ky);
			//	chunks[ky<<log_level_dx|kx].update_surfaces(kx+1>=level_dx?nullptr:&chunks[level_dx*ky+kx+1], kx-1<0?nullptr:&chunks[level_dx*ky+kx-1], ky+1>=level_dy?nullptr:&chunks[level_dx*(ky+1)+kx], ky-1<0?nullptr:&chunks[level_dx*(ky-1)+kx]);
		prof_add("Update surfaces");
#ifdef LIGHT_SYSTEM
		if(lighting)
		{
			calculate_lights_v3(true);
		//	calculate_lights_rain();
		//	calculate_lights();
			prof_add("Calculate lights");
		}
#endif
	}
	//if(dx||dy)
	//{
	//}
}
void		seed_world(unsigned long long seed)
{
	world_seed=seed;
	//srand((unsigned)world_seed);
	//unsigned caves_seed=rand()<<15|rand();
	//unsigned surface_seed=rand()<<15|rand();
	pdf_caves.seed(world_seed);
	pdf_surface.seed(~world_seed);
	xoshiro128_seed(seed, ~seed);
	//pdf_air.seed(0);//TODO: user-entered seed
	//pdf_dirt.seed(1);//
}
#ifdef LIGHT_SYSTEM
float		light_levels[16]={0};
inline float vertex_light(byte *n1, byte *n2, byte *n3, byte *n4)
{
	float sum=0;
	int idx=0;
	if(n1)
	{
		idx=maximum(*n1>>4, *n1&15);
		sum+=light_levels[idx];
	}
	if(n2)
	{
		idx=maximum(*n2>>4, *n2&15);
		sum+=light_levels[idx];
	}
	if(n3)
	{
		idx=maximum(*n3>>4, *n3&15);
		sum+=light_levels[idx];
	}
	if(n4)
	{
		idx=maximum(*n4>>4, *n4&15);
		sum+=light_levels[idx];
	}
//	float sum=(n1?light_levels[maximum(*n1>>4, *n1&15)]:0)+(n2?light_levels[maximum(*n2>>4, *n2&15)]:0)+(n3?light_levels[maximum(*n3>>4, *n3&15)]:0)+(n4?light_levels[maximum(*n4>>4, *n4&15)]:0);
	return sum*0.25f;

	//float sum=(n1?light_levels[*n1]:0)+(n2?light_levels[*n2]:0)+(n3?light_levels[*n3]:0)+(n4?light_levels[*n4]:0);
	//return sum*0.25f;

	//return light_levels[((n1?*n1:0)+(n2?*n2:0)+(n3?*n3:0)+(n4?*n4:0))>>2];//high contrast
}
//inline float vertex_light(bool f1, bool f2, bool f3, char *n1, char *n2, char *n3, char *n4, char *n5, char *n6, char *n7)
//{
//	if((int)f1|(int)f2|(int)f3)
//	{
//		float sum=(n1?light_levels[*n1]:0)+(n2?light_levels[*n2]:0)+(n3?light_levels[*n3]:0)+(n4?light_levels[*n4]:0)+(n5?light_levels[*n5]:0)+(n6?light_levels[*n6]:0)+(n7?light_levels[*n7]:0);
//		float count=float((n1!=nullptr)+(n2!=nullptr)+(n3!=nullptr)+(n4!=nullptr)+(n5!=nullptr)+(n6!=nullptr)+(n7!=nullptr));
//		return count?sum/count:0;
//	}
//	return 0;
//}
#endif

float			g_fbuf[16]={0};
namespace		GL2_2D
{
	unsigned	program=0;
	int			u_color=-1, a_vertices=-1;
	unsigned	vertex_buffer=0;
	ivec4		region, current_region;//x1, y1, x2, y2
	bool		continuous=true;
	int			pen_color=0xFF000000;
	struct		DrawInfo
	{
		int count, color;
		DrawInfo(int count, int color):count(count), color(color){}
	};
	std::vector<DrawInfo> drawInfo;
	std::vector<vec2> vertices;
	void		set_color(int color)
	{
		gl_setProgram(GL2_2D::program);
		send_color(GL2_2D::u_color, color);
	}
	void		set_region(int x1, int x2, int y1, int y2)
	{
		region.set(x1, y1, x2-x1, y2-y1);
	}
	void		use_region()
	{
		glViewport(region.x1, region.y1, region.dx, region.dy);
		current_region=region;
	}
	void		drop_region()
	{
		glViewport(0, 0, w, h);
		current_region.set(0, 0, w, h);
	}
	void		toNDC(float xs, float ys, float &xn, float &yn)
	{
		xn=(xs-current_region.x1)*2.f/current_region.dx-1;
		yn=1-(ys-current_region.y1)*2.f/current_region.dy;
	}
	//void		toNDC(int xs, int ys, float &xn, float &yn)
	//{
	//	xn=(xs-current_region.x1)*2.f/current_region.dx-1;
	//	yn=1-(ys-current_region.y1)*2.f/current_region.dy;
	//}
	void		curve_start(){GL2_2D::vertices.clear(), GL2_2D::drawInfo.clear();}
	void		curve_point(float x, float y)
	{
		vec2 ndc;
		toNDC(x, y, ndc.x, ndc.y);
		GL2_2D::vertices.push_back(ndc);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//GL2_2D::vertices.push_back(vec2(x*_2_w-1, 1-y*_2_h));
		bool increment=false;
		if(GL2_2D::drawInfo.size())
		{
			auto &di=*GL2_2D::drawInfo.rbegin();
			increment=GL2_2D::continuous&(di.color==GL2_2D::pen_color);
		}
		if(increment)
			++GL2_2D::drawInfo.rbegin()->count;
		else
			GL2_2D::drawInfo.push_back(GL2_2D::DrawInfo(1, GL2_2D::pen_color));
		GL2_2D::continuous=true;
	}
	void		draw_curve()
	{
		using namespace GL2_2D;
		gl_setProgram(program);									check(__LINE__);
		
		glEnableVertexAttribArray(a_vertices);					check(__LINE__);
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);													check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, (int)vertices.size()*sizeof(vec2), &GL2_2D::vertices[0], GL_STATIC_DRAW);	check(__LINE__);
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);											check(__LINE__);

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, &vertices[0]);								check(__LINE__);//doesn't draw

		for(int k=0, kEnd=GL2_2D::drawInfo.size(), idx=0;k<kEnd;++k)
		{
			auto &di=drawInfo[k];
			send_color(GL2_2D::u_color, di.color);								check(__LINE__);
			glDrawArrays(di.count==1?GL_POINTS:GL_LINE_STRIP, idx, di.count);	check(__LINE__);
			idx+=di.count;
		}
		glDisableVertexAttribArray(a_vertices);					check(__LINE__);
	}
	void		draw_line(float x1, float y1, float x2, float y2)//send color first
	{
		toNDC(x1, y1, g_fbuf[0], g_fbuf[1]);
		toNDC(x2, y2, g_fbuf[2], g_fbuf[3]);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//g_fbuf[0]=x1*_2_w-1, g_fbuf[1]=1-y1*_2_h;
		//g_fbuf[2]=x2*_2_w-1, g_fbuf[3]=1-y2*_2_h;
		gl_setProgram(program);					check(__LINE__);
		glEnableVertexAttribArray(a_vertices);	check(__LINE__);
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);			check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, 4<<2, g_fbuf, GL_STATIC_DRAW);	check(__LINE__);//use glBufferSubData
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);	check(__LINE__);

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	check(__LINE__);//use buffer

		glDrawArrays(GL_LINES, 0, 2);			check(__LINE__);
	}
	inline void	draw_line(int x1, int y1, int x2, int y2){draw_line((float)x1, (float)y1, (float)x2, (float)y2);}
	void		set_pixel(float x, float y)//send color first
	{
		toNDC(x, y, g_fbuf[0], g_fbuf[1]);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//g_fbuf[0]=x*_2_w-1, g_fbuf[1]=1-y*_2_h;
		gl_setProgram(program);					check(__LINE__);
		glEnableVertexAttribArray(a_vertices);	check(__LINE__);
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);			check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, 2<<2, g_fbuf, GL_STATIC_DRAW);	check(__LINE__);
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);	check(__LINE__);

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	check(__LINE__);//use buffer

		glDrawArrays(GL_POINTS, 0, 1);			check(__LINE__);
	}
	void		draw_rectangle(float x1, float x2, float y1, float y2)
	{
		float X1, X2, Y1, Y2;
		toNDC(x1, y1, X1, Y1);
		toNDC(x2, y2, X2, Y2);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//float X1=x1*_2_w-1, X2=x2*_2_w-1, Y1=1-y1*_2_h, Y2=1-y2*_2_h;
		g_fbuf[0]=X1, g_fbuf[1]=Y2;
		g_fbuf[2]=X2, g_fbuf[3]=Y2;
		g_fbuf[4]=X2, g_fbuf[5]=Y1;
		g_fbuf[6]=X1, g_fbuf[7]=Y1;
		g_fbuf[8]=X1, g_fbuf[9]=Y2;
		gl_setProgram(program);					check(__LINE__);
		
		glBindBuffer(GL_ARRAY_BUFFER, GL2_2D::vertex_buffer);			check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, 10<<2, g_fbuf, GL_STATIC_DRAW);	check(__LINE__);
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, 0);	check(__LINE__);

	//	glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	check(__LINE__);//use buffer

		glEnableVertexAttribArray(a_vertices);	check(__LINE__);
		glDrawArrays(GL_TRIANGLE_FAN, 0, 5);	check(__LINE__);
	}
	inline void	draw_rectangle(int x1, int x2, int y1, int y2){draw_rectangle((float)x1, (float)x2, (float)y1, (float)y2);}
	void		draw_rectangle_hollow(float x1, float x2, float y1, float y2)
	{
		float X1, X2, Y1, Y2;
		toNDC(x1, y1, X1, Y1);
		toNDC(x2, y2, X2, Y2);
		//float _2_w=2.f/w, _2_h=2.f/h;
		//float X1=x1*_2_w-1, X2=x2*_2_w-1, Y1=1-y1*_2_h, Y2=1-y2*_2_h;
		g_fbuf[0]=X1, g_fbuf[1]=Y1;
		g_fbuf[2]=X2, g_fbuf[3]=Y1;
		g_fbuf[4]=X2, g_fbuf[5]=Y2;
		g_fbuf[6]=X1, g_fbuf[7]=Y2;
		g_fbuf[8]=X1, g_fbuf[9]=Y1;

		glBindBuffer(GL_ARRAY_BUFFER, 0);										check(__LINE__);
		glVertexAttribPointer(a_vertices, 2, GL_FLOAT, GL_FALSE, 0, g_fbuf);	check(__LINE__);

		glEnableVertexAttribArray(a_vertices);	check(__LINE__);
		glDrawArrays(GL_LINE_STRIP, 0, 5);		check(__LINE__);
	}
	void		draw_rectangle_hollow(int x1, int x2, int y1, int y2){draw_rectangle_hollow((float)x1, (float)x2, (float)y1, (float)y2);}
}
namespace		GL2_L3D
{
	unsigned	program=0;
	int			a_vertices=-1, a_normals=-1, u_vpmatrix=-1, u_modelmatrix=-1, u_normalmatrix=-1, u_objectcolor=-1, u_lightcolor=-1, u_lightpos=-1, u_viewpos=-1;
}
unsigned		sphereVBO=0, sphereEBO=0;

typedef std::vector<vec3> TxVertices;//VVVTTL: vec3 vertex, vec2 texcoord, float light
namespace		GL2_BB3D
{
	unsigned	program=0;
	int			a_vertices=-1, a_texcoord=-1, a_light=-1, u_matrix=-1, u_texture=-1;
//	unsigned	buffer=0;
	Chunk		*ch;
//	std::vector<BBQuadInfo> bbDrawInfo;
	std::vector<TxVertices> bbVertices;
	int			transparent_start_idx=-1;
	void		begin_chunk(Chunk *ch)
	{
		GL2_BB3D::ch=ch;
		ch->glDrawInfo.clear();
	//	GL2_BB3D::bbDrawInfo.clear();
		ch->glDrawInfo.reserve(nTextures);
		bbVertices.reserve(nTextures);
		for(int k=0;k<nTextures;++k)
			ch->glDrawInfo.push_back(BBQuadInfo(textures[k+1])), bbVertices.push_back(TxVertices());
	}
//	void		begin_transparent(){GL2_BB3D::transparent_start_idx=GL2_BB3D::bbDrawInfo.size();}//must be called once
	void		push_face_initiate(TxVertices *&tv, BBQuadInfo *&bbdi, unsigned tx_id)
	{
		for(int k=0, kEnd=GL2_BB3D::ch->glDrawInfo.size();k<kEnd;++k)//find drawinfo structure with this texture id
		{
			auto &bbdi2=ch->glDrawInfo[k];
			if(bbdi2.tx_id==tx_id)
			{
				bbdi=&bbdi2, tv=&bbVertices[k];
				break;
			}
		}
	}
	void		push_face_down(int x, int y, int z, int dx, int dy, unsigned tx_id, vec4 const &light, bool flipped)//brick blaster
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		push_face_initiate(tv, bbdi, tx_id);
		if(bbdi&&tv)
		{
			if(flipped)
			{
				tv->push_back(vec3((float)(x+dx),	(float)(y+dy),	(float)z)), tv->push_back(vec3(1, 1, light.z));//DNE	//down face
				tv->push_back(vec3((float)(x+dx),	(float)y,		(float)z)), tv->push_back(vec3(1, 0, light.y));//DSE
				tv->push_back(vec3((float)x,		(float)y,		(float)z)), tv->push_back(vec3(0, 0, light.x));//DSW
				tv->push_back(vec3((float)x,		(float)(y+dy),	(float)z)), tv->push_back(vec3(0, 1, light.w));//DNW
			}
			else
			{
				tv->push_back(vec3((float)x,		(float)(y+dy),	(float)z)), tv->push_back(vec3(0, 1, light.w));//DNW	//down face
				tv->push_back(vec3((float)(x+dx),	(float)(y+dy),	(float)z)), tv->push_back(vec3(1, 1, light.z));//DNE
				tv->push_back(vec3((float)(x+dx),	(float)y,		(float)z)), tv->push_back(vec3(1, 0, light.y));//DSE
				tv->push_back(vec3((float)x,		(float)y,		(float)z)), tv->push_back(vec3(0, 0, light.x));//DSW
			}
			bbdi->count+=4;
		}
	}
	void		push_face_up(int x, int y, int z, int dx, int dy, unsigned tx_id, vec4 const &light, bool flipped)
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		push_face_initiate(tv, bbdi, tx_id);
		if(bbdi&&tv)
		{
			if(flipped)
			{
				tv->push_back(vec3((float)(x+dx),	(float)y,		(float)z)), tv->push_back(vec3(1, 0, light.y));//USE	//up face
				tv->push_back(vec3((float)(x+dx),	(float)(y+dy),	(float)z)), tv->push_back(vec3(1, 1, light.z));//UNE
				tv->push_back(vec3((float)x,		(float)(y+dy),	(float)z)), tv->push_back(vec3(0, 1, light.w));//UNW
				tv->push_back(vec3((float)x,		(float)y,		(float)z)), tv->push_back(vec3(0, 0, light.x));//USW
			}
			else
			{
				tv->push_back(vec3((float)x,		(float)y,		(float)z)), tv->push_back(vec3(0, 0, light.x));//USW	//up face
				tv->push_back(vec3((float)(x+dx),	(float)y,		(float)z)), tv->push_back(vec3(1, 0, light.y));//USE
				tv->push_back(vec3((float)(x+dx),	(float)(y+dy),	(float)z)), tv->push_back(vec3(1, 1, light.z));//UNE
				tv->push_back(vec3((float)x,		(float)(y+dy),	(float)z)), tv->push_back(vec3(0, 1, light.w));//UNW
			}
			bbdi->count+=4;
		}
	}
	void		push_face_north(int x, int y, int z, int dx, int dz, unsigned tx_id, vec4 const &light, bool flipped)
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		push_face_initiate(tv, bbdi, tx_id);
		if(bbdi&&tv)
		{
			if(flipped)
			{
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.z));//UNE	//north face
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.y));//DNE
				tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.x));//DNW
				tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.w));//UNW
			}
			else
			{
				tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.w));//UNW	//north face
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.z));//UNE
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.y));//DNE
				tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.x));//DNW
			}
			bbdi->count+=4;
		}
	}
	void		push_face_south(int x, int y, int z, int dx, int dz, unsigned tx_id, vec4 const &light, bool flipped)
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		push_face_initiate(tv, bbdi, tx_id);
		if(bbdi&&tv)
		{
			if(flipped)
			{
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DSE	//south face
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//USE
				tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USW
				tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSW
			}
			else
			{
				tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSW	//south face
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DSE
				tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//USE
				tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USW
			}
			bbdi->count+=4;
		}
	}
	void		push_face_east(int x, int y, int z, int dy, int dz, unsigned tx_id, vec4 const &light, bool flipped)
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		push_face_initiate(tv, bbdi, tx_id);
		if(bbdi&&tv)
		{
			if(flipped)
			{
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DNE	//east face
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//UNE
				tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USE
				tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSE
			}
			else
			{
				tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSE	//east face
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DNE
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//UNE
				tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USE
			}
			bbdi->count+=4;
		}
	}
	void		push_face_west(int x, int y, int z, int dy, int dz, unsigned tx_id, vec4 const &light, bool flipped)
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		push_face_initiate(tv, bbdi, tx_id);
		if(bbdi&&tv)
		{
			if(flipped)
			{
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.z));//UNW	//west face
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(0, 0, light.y));//DNW
				tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(1, 0, light.x));//DSW
				tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.w));//USW
			}
			else
			{
				tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.w));//USW	//west face
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.z));//UNW
				tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(0, 0, light.y));//DNW
				tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(1, 0, light.x));//DSW
			}
			bbdi->count+=4;
		}
	}
	void		push_quad(int x, int y, int z, int dx, int dy, int dz, unsigned tx_id, bool reversed, vec4 const &light)//brick blaster
	{
		TxVertices *tv=nullptr;
		BBQuadInfo *bbdi=nullptr;
		for(int k=0, kEnd=ch->glDrawInfo.size();k<kEnd;++k)//find drawinfo structure with this texture id
		{
			auto &bbdi2=ch->glDrawInfo[k];
			if(bbdi2.tx_id==tx_id)
			{
				bbdi=&bbdi2, tv=&bbVertices[k];
				break;
			}
		}
		if(bbdi&&tv)
		{
			if(!dz)
			{
				if(reversed)
				{
					tv->push_back(vec3((float)x,		(float)(y+dy),	(float)z)), tv->push_back(vec3(0, 1, light.w));//DNW	//down face
					tv->push_back(vec3((float)(x+dx),	(float)(y+dy),	(float)z)), tv->push_back(vec3(1, 1, light.z));//DNE
					tv->push_back(vec3((float)(x+dx),	(float)y,		(float)z)), tv->push_back(vec3(1, 0, light.y));//DSE
					tv->push_back(vec3((float)x,		(float)y,		(float)z)), tv->push_back(vec3(0, 0, light.x));//DSW
				}
				else
				{
					tv->push_back(vec3((float)x,		(float)y,		(float)z)), tv->push_back(vec3(0, 0, light.x));//USW	//up face
					tv->push_back(vec3((float)(x+dx),	(float)y,		(float)z)), tv->push_back(vec3(1, 0, light.y));//USE
					tv->push_back(vec3((float)(x+dx),	(float)(y+dy),	(float)z)), tv->push_back(vec3(1, 1, light.z));//UNE
					tv->push_back(vec3((float)x,		(float)(y+dy),	(float)z)), tv->push_back(vec3(0, 1, light.w));//UNW
				}
			}
			else if(!dy)
			{
				if(reversed)
				{
					tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.w));//UNW	//north face
					tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.z));//UNE
					tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.y));//DNE
					tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.x));//DNW
					//tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//UNW	//north face
					//tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//UNE
					//tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DNE
					//tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.x));//DNW
				}
				else
				{
					tv->push_back(vec3((float)x,		(float)y,	(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSW	//south face
					tv->push_back(vec3((float)(x+dx),	(float)y,	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DSE
					tv->push_back(vec3((float)(x+dx),	(float)y,	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//USE
					tv->push_back(vec3((float)x,		(float)y,	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USW
				}
			}
			else//dx==0
			{
				if(reversed)
				{
					tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.w));//USW	//west face
					tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.z));//UNW
					tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(0, 0, light.y));//DNW
					tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(1, 0, light.x));//DSW
					//tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USW	//west face
					//tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//UNW
					//tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DNW
					//tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSW
				}
				else
				{
					tv->push_back(vec3((float)x,	(float)y,		(float)z		)), tv->push_back(vec3(0, 0, light.x));//DSE	//east face
					tv->push_back(vec3((float)x,	(float)(y+dy),	(float)z		)), tv->push_back(vec3(1, 0, light.y));//DNE
					tv->push_back(vec3((float)x,	(float)(y+dy),	(float)(z+dz)	)), tv->push_back(vec3(1, 1, light.z));//UNE
					tv->push_back(vec3((float)x,	(float)y,		(float)(z+dz)	)), tv->push_back(vec3(0, 1, light.w));//USE
				}
			}
			bbdi->count+=4;
		//	drawInfo.push_back(VertexInfo(4, tx_id));
		}
	}
	void		send_chunk()//send chunk vertices to GPU
	{
		int totalsize=0;
		for(int k=0, kEnd=GL2_BB3D::bbVertices.size();k<kEnd;++k)
			totalsize+=GL2_BB3D::bbVertices[k].size();
		std::vector<vec3> vertices2;
		vertices2.reserve(totalsize);
		for(int k=GL2_BB3D::bbVertices.size()-1;k>=0;--k)
	//	for(int k=0, kEnd=GL2_BB3D::bbVertices.size();k<kEnd;++k)
		{
			//if(k==nTextures-nTransparentTextures)
			//	ch->transparent_start_idx=
			auto &v3=bbVertices[k];
			vertices2.insert(vertices2.end(), v3.begin(), v3.end());
		}
		glBindBuffer(GL_ARRAY_BUFFER, ch->gl_buffer);														check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, (int)vertices2.size()*sizeof(vec3), &vertices2[0], GL_STATIC_DRAW);	check(__LINE__);
		bbVertices.clear();
	}
	void		draw(mat4 const &mVP)//draws all chunks
	{
		gl_setProgram(GL2_BB3D::program);		check(__LINE__);

		glUniformMatrix4fv(GL2_BB3D::u_matrix, 1, GL_FALSE, mVP.data());						check(__LINE__);

		std::vector<int> p_idx(chunks.size());
		for(int c=0, cEnd=chunks.size();c<cEnd;++c)//draw all opaque quads
		{
			auto &ch=chunks[c];
			glBindBuffer(GL_ARRAY_BUFFER, ch.gl_buffer);											check(__LINE__);
			glVertexAttribPointer(GL2_BB3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 6<<2, 0);			check(__LINE__);
			glVertexAttribPointer(GL2_BB3D::a_texcoord, 3, GL_FLOAT, GL_FALSE, 6<<2, (void*)(3<<2));check(__LINE__);
			glVertexAttribPointer(GL2_BB3D::a_light, 3, GL_FLOAT, GL_FALSE, 6<<2, (void*)(5<<2));	check(__LINE__);
			glEnableVertexAttribArray(GL2_BB3D::a_vertices);										check(__LINE__);
			glEnableVertexAttribArray(GL2_BB3D::a_texcoord);										check(__LINE__);
			glEnableVertexAttribArray(GL2_BB3D::a_light);											check(__LINE__);
			for(int k=ch.glDrawInfo.size()-1;k>=nTransparentTextures;--k)
		//	for(int k=0, kEnd=nTextures-nTransparentTextures;k<kEnd;++k)
			{
				auto &di=ch.glDrawInfo[k];
				if(di.count)
				{
					glBindTexture(GL_TEXTURE_2D, di.tx_id);			check(__LINE__);
					glDrawArrays(GL_QUADS, p_idx[c], di.count);		check(__LINE__);
					p_idx[c]+=di.count;
				}
			}
		}
		glDepthMask(GL_FALSE);
		for(int c=0, cEnd=chunks.size();c<cEnd;++c)//draw all transparent quads
		{
			auto &ch=chunks[c];
			glBindBuffer(GL_ARRAY_BUFFER, ch.gl_buffer);											check(__LINE__);
			glVertexAttribPointer(GL2_BB3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 6<<2, 0);			check(__LINE__);
			glVertexAttribPointer(GL2_BB3D::a_texcoord, 3, GL_FLOAT, GL_FALSE, 6<<2, (void*)(3<<2));check(__LINE__);
			glVertexAttribPointer(GL2_BB3D::a_light, 3, GL_FLOAT, GL_FALSE, 6<<2, (void*)(5<<2));	check(__LINE__);
			glEnableVertexAttribArray(GL2_BB3D::a_vertices);										check(__LINE__);
			glEnableVertexAttribArray(GL2_BB3D::a_texcoord);										check(__LINE__);
			glEnableVertexAttribArray(GL2_BB3D::a_light);											check(__LINE__);
			for(int k=nTransparentTextures-1;k>=0;--k)
		//	for(int k=nTextures-nTransparentTextures, kEnd=ch.glDrawInfo.size();k<kEnd;++k)
			{
				auto &di=ch.glDrawInfo[k];
				if(di.count)
				{
					glBindTexture(GL_TEXTURE_2D, di.tx_id);			check(__LINE__);
					glDrawArrays(GL_QUADS, p_idx[c], di.count);		check(__LINE__);
					p_idx[c]+=di.count;
				}
			}
		}
		glDepthMask(GL_TRUE);
		glDisableVertexAttribArray(GL2_BB3D::a_vertices);	check(__LINE__);
		glDisableVertexAttribArray(GL2_BB3D::a_texcoord);	check(__LINE__);
		glDisableVertexAttribArray(GL2_BB3D::a_light);		check(__LINE__);
	//	glBindBuffer(GL_ARRAY_BUFFER, 0);					check(__LINE__);//for immediate mode
	}
}
struct			VertexInfo//for 3D program
{
	int count;
	unsigned type;//GL_POINTS/GL_LINES/GL_LINE_STRIP/GL_TRIANGLES/GL_TRIANGLE_FAN/GL_QUADS/...
	char textured;//0: color, 1: tx pointer, 2: use color as texture id
	int color, *tx, txw, txh;
	VertexInfo(int count, unsigned tx_id):count(count), type(GL_QUADS), textured(2), color(tx_id), tx(nullptr), txw(), txh(){}//brick blaster vertex
	VertexInfo(int count, unsigned type, int color):count(count), type(type), textured(false), color(color), tx(nullptr), txw(0), txh(0){}
	VertexInfo(int count, unsigned type, int *tx, int txw, int txh):count(count), type(type), textured(true), color(0), tx(tx), txw(txw), txh(txh){}
};
namespace		GL2_3D
{
	unsigned	program=0;
	int			a_vertices=-1, a_texcoord=-1, u_matrix=-1, u_pointSize=-1, u_texture=-1, u_isTextured=-1, u_color=-1;
	unsigned	vertex_buffer=0, texcoord_buffer=0;
	int			transparent_start_idx=-1, bb_start_idx=-1, bb_tstart_idx=-1;
	std::vector<vec3> vertices;
	std::vector<vec2> texcoord;
	std::vector<VertexInfo> drawInfo;
	void		begin(){drawInfo.clear();}
	void		begin_transparent(){GL2_3D::transparent_start_idx=GL2_3D::drawInfo.size();}
	void		push_square(float x1, float x2, float y1, float y2, float z, int *tx, int txw, int txh)
	{
		vertices.push_back(vec3(x1, y1, z)), texcoord.push_back(vec2(0, 0));
		vertices.push_back(vec3(x2, y1, z)), texcoord.push_back(vec2(1, 0));
		vertices.push_back(vec3(x2, y2, z)), texcoord.push_back(vec2(1, 1));
		vertices.push_back(vec3(x1, y2, z)), texcoord.push_back(vec2(0, 1));
		drawInfo.push_back(VertexInfo(4, GL_QUADS, tx, txw, txh));
	}
	void		push_triangle(vec3 const &a, vec3 const &b, vec3 const &c, int color)
	{
		vertices.push_back(a), texcoord.push_back(vec2(0, 0));
		vertices.push_back(b), texcoord.push_back(vec2(0, 0));
		vertices.push_back(c), texcoord.push_back(vec2(0, 0));
		bool increment=false;
		if(drawInfo.size())
		{
			auto &di=*drawInfo.rbegin();
			increment=di.type==GL_TRIANGLES&&di.color==color;
		}
		if(increment)
			drawInfo.rbegin()->count+=3;
		else
			drawInfo.push_back(VertexInfo(3, GL_TRIANGLES, color));
	}
	void		push_line_segment(vec3 const &p1, vec3 const &p2, int color)
	{
		vertices.push_back(p1), texcoord.push_back(vec2(0, 0));
		vertices.push_back(p2), texcoord.push_back(vec2(0, 0));
		bool increment=false;
		if(drawInfo.size())
		{
			auto &di=*drawInfo.rbegin();
			increment=di.type==GL_LINES&&di.color==color;
		}
		if(increment)
			drawInfo.rbegin()->count+=2;
		else
			drawInfo.push_back(VertexInfo(2, GL_LINES, color));
	}
	void push_point(vec3 const &p, int color)
	{
		vertices.push_back(p), texcoord.push_back(vec2(0, 0));
		bool increment=false;
		if(drawInfo.size())
		{
			auto &di=*drawInfo.rbegin();
			increment=di.type==GL_POINTS&&di.color==color;
		}
		if(increment)
			drawInfo.rbegin()->count+=1;
		else
			drawInfo.push_back(VertexInfo(1, GL_POINTS, color));
	}
	void end()
	{
		glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::vertex_buffer);													check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, (int)vertices.size()*sizeof(vec3), &GL2_3D::vertices[0], GL_STATIC_DRAW);	check(__LINE__);
		glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::texcoord_buffer);													check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, (int)texcoord.size()*sizeof(vec2), &GL2_3D::texcoord[0], GL_STATIC_DRAW);	check(__LINE__);
		vertices.clear(), texcoord.clear();
	}
	void draw(mat4 const &mVP)
	{
		gl_setProgram(GL2_3D::program);		check(__LINE__);

		glUniformMatrix4fv(GL2_3D::u_matrix, 1, GL_FALSE, mVP.data());			check(__LINE__);

		static unsigned tx_id=0;
		if(!tx_id)
			{glGenTextures(1, &tx_id);			check(__LINE__);}

		glEnableVertexAttribArray(GL2_3D::a_vertices);							check(__LINE__);
		glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);							check(__LINE__);
		glVertexAttribPointer(GL2_3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 0, 0);	check(__LINE__);
		glEnableVertexAttribArray(GL2_3D::a_texcoord);							check(__LINE__);
		glBindBuffer(GL_ARRAY_BUFFER, GL2_3D::texcoord_buffer);					check(__LINE__);
		glVertexAttribPointer(GL2_3D::a_texcoord, 2, GL_FLOAT, GL_FALSE, 0, 0);	check(__LINE__);
		for(int k=0, p_idx=0, kEnd=drawInfo.size();k<kEnd;++k)
		{
			auto &di=GL2_3D::drawInfo[k];
			if(k==transparent_start_idx)
				glDepthMask(GL_FALSE);
			if(k==bb_start_idx)
				glDepthMask(GL_TRUE);
			if(k==bb_tstart_idx)
				glDepthMask(GL_FALSE);
			glUniform1i(GL2_3D::u_isTextured, di.textured);
			send_color(u_color, di.textured?0xFFFF00FF:di.color);				check(__LINE__);//dummy color if textured
			if(di.textured==1)
			{
				glBindTexture(GL_TEXTURE_2D, tx_id);															check(__LINE__);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);								check(__LINE__);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);								check(__LINE__);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, di.txw, di.txh, 0, GL_RGBA, GL_UNSIGNED_BYTE, di.tx);	check(__LINE__);//TODO: send textures to GPU with the vertices
				select_texture(tx_id, GL2_3D::u_texture);														check(__LINE__);
			}
			else if(di.textured==2)
			{
				glBindTexture(GL_TEXTURE_2D, di.color);			check(__LINE__);
			}
			else
				{select_texture(tx_id, GL2_3D::u_texture);		check(__LINE__);}//dummy texture if not textured
			glDrawArrays(di.type, p_idx, di.count);				check(__LINE__);
			p_idx+=di.count;
		}
		glDisableVertexAttribArray(GL2_3D::a_vertices);	check(__LINE__);
		glDisableVertexAttribArray(GL2_3D::a_texcoord);	check(__LINE__);
	//	glBindBuffer(GL_ARRAY_BUFFER, 0);				check(__LINE__);
		glDepthMask(GL_TRUE);
	}
	unsigned dummy_tx_id=0;
	void draw_line(mat4 const &mVP, vec3 const &p1, vec3 const &p2, int color)//TODO: pre-calculate matrices once per frame
	{
		gl_setProgram(GL2_3D::program);		check(__LINE__);

		glUniformMatrix4fv(GL2_3D::u_matrix, 1, GL_FALSE, mVP.data());					check(__LINE__);

		if(!dummy_tx_id)
		{
			int txw=1, txh=1, tx=0;//dummy 1x1 texture
			glGenTextures(1, &dummy_tx_id);															check(__LINE__);
			glBindTexture(GL_TEXTURE_2D, dummy_tx_id);												check(__LINE__);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);						check(__LINE__);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);						check(__LINE__);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txw, txh, 0, GL_RGBA, GL_UNSIGNED_BYTE, &tx);	check(__LINE__);//TODO: send textures to GPU with the vertices
		}
		
		float vertices[]=
		{
			p1.x, p1.y, p1.z,
			p2.x, p2.y, p2.z,
		};
		const float texcoord[4]={0};
		glBindBuffer(GL_ARRAY_BUFFER, 0);												check(__LINE__);
		glVertexAttribPointer(GL2_3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 0, vertices);	check(__LINE__);
		glVertexAttribPointer(GL2_3D::a_texcoord, 2, GL_FLOAT, GL_FALSE, 0, texcoord);	check(__LINE__);
		glEnableVertexAttribArray(GL2_3D::a_vertices);									check(__LINE__);
		glEnableVertexAttribArray(GL2_3D::a_texcoord);									check(__LINE__);

		glUniform1i(GL2_3D::u_isTextured, 0);			check(__LINE__);
		send_color(u_color, color);						check(__LINE__);
		select_texture(dummy_tx_id, GL2_3D::u_texture);	check(__LINE__);//dummy texture
		glDrawArrays(GL_LINES, 0, 2);					check(__LINE__);

		glDisableVertexAttribArray(GL2_3D::a_vertices);	check(__LINE__);
		glDisableVertexAttribArray(GL2_3D::a_texcoord);	check(__LINE__);
	}
	void draw_AABB(mat4 const &mVP, float x1, float x2, float y1, float y2, float z1, float z2, int color)
	{
		gl_setProgram(GL2_3D::program);		check(__LINE__);

		glUniformMatrix4fv(GL2_3D::u_matrix, 1, GL_FALSE, mVP.data());					check(__LINE__);

		if(!dummy_tx_id)
		{
			int txw=1, txh=1, tx=0;//dummy 1x1 texture
			glGenTextures(1, &dummy_tx_id);															check(__LINE__);
			glBindTexture(GL_TEXTURE_2D, dummy_tx_id);												check(__LINE__);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);						check(__LINE__);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);						check(__LINE__);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txw, txh, 0, GL_RGBA, GL_UNSIGNED_BYTE, &tx);	check(__LINE__);//TODO: send textures to GPU with the vertices
		}
		
		float vertices[]=//24 vertices
		{
			x1, y1, z1,		x2, y1, z1,
			x2, y1, z1,		x2, y2, z1,
			x2, y2, z1,		x1, y2, z1,
			x1, y2, z1,		x1, y1, z1,

			x1, y1, z1,		x1, y1, z2,
			x2, y1, z1,		x2, y1, z2,
			x2, y2, z1,		x2, y2, z2,
			x1, y2, z1,		x1, y2, z2,

			x1, y1, z2,		x2, y1, z2,
			x2, y1, z2,		x2, y2, z2,
			x2, y2, z2,		x1, y2, z2,
			x1, y2, z2,		x1, y1, z2,
		};
		const float texcoord[24<<1]={0};
		glBindBuffer(GL_ARRAY_BUFFER, 0);												check(__LINE__);
		glVertexAttribPointer(GL2_3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 0, vertices);	check(__LINE__);
		glVertexAttribPointer(GL2_3D::a_texcoord, 2, GL_FLOAT, GL_FALSE, 0, texcoord);	check(__LINE__);
		glEnableVertexAttribArray(GL2_3D::a_vertices);									check(__LINE__);
		glEnableVertexAttribArray(GL2_3D::a_texcoord);									check(__LINE__);

		glUniform1i(GL2_3D::u_isTextured, 0);			check(__LINE__);
		send_color(u_color, color);						check(__LINE__);
		select_texture(dummy_tx_id, GL2_3D::u_texture);	check(__LINE__);//dummy texture
		glDrawArrays(GL_LINES, 0, 24);					check(__LINE__);

		glDisableVertexAttribArray(GL2_3D::a_vertices);	check(__LINE__);
		glDisableVertexAttribArray(GL2_3D::a_texcoord);	check(__LINE__);
	}
}
unsigned		programID=0;
int				attribute_coord3d=0, attribute_texcoord=0, uniform_mvp=0, uniform_texture=0, uniform_fade=0;
//TODO: namespace Resources
#if 1
const char sf10_height=16, sf8_height=12, sf7_height=10, sf6_height=8;
const char sf10_widths[]=
{
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,4,4,6,8,8,11,9,4,
	4,4,6,8,4,4,4,4,8,8,
	8,8,8,8,8,8,8,8,4,4,
	8,8,8,8,14,8,10,9,10,9,
	8,10,10,4,7,9,8,12,10,10,
	9,10,10,9,8,10,8,14,9,10,
	9,4,4,4,5,8,5,8,8,7,
	8,8,4,8,8,4,4,7,4,12,
	8,8,8,8,5,8,4,8,8,10,
	8,8,8,5,4,5,5,0
};
const char sf8_widths[]=
{
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,2,3,4,6,6,9,7,2,
	3,3,4,6,3,5,3,3,6,4,
	6,6,6,6,6,6,6,6,3,3,
	6,6,6,6,10,7,7,7,7,6,
	6,8,7,2,5,7,6,8,7,8,
	7,8,7,7,6,7,7,9,7,7,
	6,3,3,3,3,6,3,6,6,6,
	6,6,3,6,6,2,2,5,2,8,
	6,6,6,6,3,5,3,6,5,7,
	5,5,5,3,3,3,3,0,
};
const char sf7_widths[]=
{
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,2,2,4,6,4,5,5,2,
	3,3,3,4,2,3,2,3,5,3,
	5,5,5,4,5,4,5,5,2,2,
	4,4,4,5,8,6,6,6,6,5,
	5,6,5,2,4,5,4,8,6,6,
	5,6,6,5,5,6,6,8,6,6,
	5,3,3,3,4,5,3,4,5,4,
	5,4,3,5,4,2,2,4,2,6,
	4,5,5,5,3,4,3,4,4,6,
	4,4,4,3,2,3,4,0,
};
char sf6_widths[]=
{
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,1,2,3,5,4,5,5,2,
	3,3,3,4,2,3,2,3,4,3,
	4,4,4,4,4,4,4,4,2,2,
	3,3,3,3,7,5,5,5,5,4,
	4,5,5,2,4,5,4,6,5,5,
	5,5,5,5,4,5,4,6,4,4,
	4,3,3,3,4,4,3,4,4,4,
	4,4,3,4,4,2,2,4,2,6,
	4,4,4,4,3,3,3,4,4,5,
	4,4,4,3,2,3,4,0,
};
int txw=128, txh=256;
namespace		resources
{
	const unsigned char sf10[]=//system font size 10+
	{
		25,209,255,159,199,24,205,243,60,99,80,217,216,252,217,216,108,108,254,108,108,50,182,140,55,207,131,7,195,179,199,96,86,29,188,177,205,108,195,13,176,195,54,179,141,61,184,204,41,199,
		102,179,113,28,219,249,204,159,17,250,207,76,123,219,182,109,155,49,211,179,109,219,182,61,140,21,243,51,158,148,178,99,24,230,103,24,6,59,99,15,72,241,89,145,199,16,51,51,155,153,
		205,204,100,76,189,121,158,231,121,158,103,47,67,202,207,204,204,204,156,49,245,230,25,198,24,99,24,254,140,169,55,195,48,7,195,240,236,101,76,97,156,231,109,123,254,97,56,99,250,15,
		195,240,205,48,60,123,25,83,111,30,134,111,158,231,217,203,152,254,48,134,49,12,99,24,70,198,212,155,231,217,155,231,121,246,50,166,222,60,207,179,15,195,179,151,18,124,128,71,77,218,
		0,96,47,100,9,99,140,49,24,12,6,199,140,126,0,248,67,150,131,193,96,48,198,24,67,198,212,155,103,24,99,24,128,97,34,220,240,192,57,6,102,110,179,61,219,155,189,205,219,108,
		123,6,192,113,240,129,65,197,192,224,225,33,49,51,243,27,30,206,160,254,13,15,15,255,13,15,15,15,255,101,78,121,230,225,96,48,24,12,205,188,12,234,207,216,240,240,240,240,240,240,
		216,79,230,244,63,24,12,126,131,193,96,240,207,156,254,7,131,193,111,48,24,12,6,25,84,62,227,193,129,129,249,225,97,99,222,12,234,240,240,240,240,255,240,240,240,240,112,70,244,255,
		255,49,166,48,12,195,48,12,207,179,151,65,29,155,217,120,56,120,216,152,25,27,206,156,14,6,131,193,96,48,24,12,254,25,213,129,7,62,252,240,231,159,223,123,239,153,103,206,160,14,
		31,63,63,111,111,207,207,143,15,103,80,121,204,134,135,135,135,135,135,205,120,50,168,127,195,195,195,195,127,3,3,3,3,25,84,30,179,225,225,225,225,225,121,51,254,76,234,159,97,195,
		134,13,251,51,108,216,176,193,153,83,223,120,60,112,224,192,227,177,15,131,250,143,129,129,129,129,129,129,129,129,145,65,29,30,30,30,30,30,30,30,54,227,193,160,14,15,155,153,153,145,
		240,240,96,96,96,92,135,225,97,120,24,54,207,204,51,150,134,231,193,48,48,12,12,3,147,58,120,176,49,54,56,112,176,49,54,120,48,70,117,224,129,13,99,6,15,24,96,128,1,6,
		24,48,169,255,1,3,3,3,3,3,3,3,3,254,103,166,191,109,219,182,237,24,226,153,25,51,51,102,102,204,244,109,219,182,109,31,83,67,220,198,99,252,71,104,206,56,102,231,205,188,
		61,207,126,198,116,24,134,111,158,231,121,254,98,118,222,60,12,195,236,101,76,97,24,246,231,121,158,103,63,102,231,205,255,48,204,30,134,148,205,222,204,204,44,102,234,207,243,60,207,62,
		60,123,25,211,97,24,190,121,158,231,121,206,136,62,255,63,102,218,134,109,219,182,151,49,29,134,97,222,158,227,217,230,140,232,255,255,199,244,254,155,121,230,153,103,158,121,230,152,221,55,
		207,243,60,207,49,59,111,158,231,121,246,98,166,223,60,207,243,252,13,195,16,51,245,231,121,158,103,31,134,225,24,221,223,204,204,196,236,188,121,240,224,217,3,81,154,189,153,153,113,204,
		110,158,231,121,158,125,24,222,240,176,153,25,15,6,6,76,111,224,153,109,179,205,63,102,152,129,225,13,155,241,96,240,152,13,195,80,135,135,205,204,120,120,48,48,24,12,98,118,63,140,
		49,198,240,103,168,217,204,204,198,204,204,56,35,253,255,255,63,134,122,204,204,140,205,204,108,48,37,183,3,
	};
	const unsigned char sf8[]=//small fonts size 8
	{
		145,184,111,100,164,141,212,41,245,213,87,42,82,72,92,53,172,58,146,185,27,142,32,8,130,56,28,49,52,146,196,72,163,200,141,132,28,145,180,170,50,34,153,170,90,68,228,67,85,132,
		124,66,146,34,22,163,248,100,196,136,160,106,85,164,208,197,24,99,140,46,50,200,147,36,25,41,116,17,34,34,194,143,20,186,8,25,132,209,69,10,17,83,41,125,132,34,133,63,132,7,
		97,116,145,66,23,195,139,49,186,72,225,135,8,17,66,72,164,208,197,232,98,140,46,82,232,98,244,16,70,151,138,138,169,176,56,67,135,36,33,132,41,154,135,207,208,133,16,146,68,164,
		208,69,136,8,1,36,131,200,35,148,89,89,153,38,192,19,49,100,24,73,146,254,48,140,24,126,97,248,133,97,248,69,12,189,48,8,130,32,244,34,134,95,24,134,97,24,126,145,194,31,
		194,11,33,252,72,225,15,225,133,16,66,228,208,23,12,4,226,193,160,47,98,24,134,225,31,134,97,24,9,252,143,16,34,34,98,166,69,12,195,40,57,142,36,10,35,133,33,132,16,66,
		248,145,195,241,184,90,77,38,131,193,136,225,56,150,101,154,198,113,228,208,23,12,6,131,193,160,47,98,248,133,97,248,5,65,16,185,244,5,131,193,96,52,236,3,70,12,191,48,12,191,
		48,12,35,134,94,24,120,32,24,122,145,194,79,8,33,132,144,136,97,24,134,97,24,134,94,196,48,12,67,73,146,24,38,130,24,24,24,40,164,165,69,66,18,49,12,67,137,97,164,48,
		132,28,6,131,34,10,2,129,64,68,10,63,68,68,132,240,35,146,87,85,103,228,170,213,136,228,170,234,144,17,21,106,241,71,68,50,85,229,160,143,62,82,24,66,120,49,198,151,170,114,
		49,68,23,41,132,16,250,24,163,79,85,185,248,131,139,8,90,87,165,234,124,140,209,195,23,41,12,33,180,25,99,204,196,245,51,145,253,35,132,17,145,53,149,145,192,255,212,213,47,153,
		76,38,83,85,47,198,24,83,85,46,198,232,82,117,47,198,248,66,72,213,249,24,163,135,48,37,117,85,138,202,195,240,66,98,93,153,170,138,49,102,155,162,202,76,179,84,85,172,85,169,
		20,85,166,101,166,232,50,83,73,34,69,245,36,241,17,73,213,84,35,145,255,35,146,85,86,65,70,60,
	};
	const unsigned char sf7[]=//small fonts size 7
	{
		145,176,55,50,210,70,202,212,87,234,171,200,156,62,249,138,144,205,146,52,71,200,164,164,75,35,33,71,4,173,42,35,130,169,106,17,145,15,205,232,202,133,156,146,56,23,49,34,166,86,
		17,50,203,204,180,136,152,171,70,200,44,209,226,35,100,150,132,105,17,50,50,171,167,200,216,51,114,17,50,139,203,180,200,216,41,73,17,50,75,203,180,8,153,165,99,90,40,40,67,81,
		57,51,69,69,12,205,28,207,76,69,84,34,100,150,36,32,145,59,142,40,89,205,9,56,145,50,66,148,139,49,82,246,226,139,241,69,202,92,12,33,186,72,217,139,49,198,23,33,251,184,
		136,143,144,125,92,68,68,202,92,12,57,250,8,89,230,103,102,36,236,143,140,145,100,23,33,203,154,169,140,140,37,73,30,57,27,143,171,213,100,50,82,54,103,173,57,71,202,92,140,49,
		186,8,217,101,94,68,164,204,197,88,147,141,148,189,24,95,140,17,50,75,161,52,72,217,39,132,16,18,41,139,49,198,232,34,101,49,42,69,72,228,44,24,84,169,40,148,72,89,84,132,
		168,24,41,139,138,16,66,34,100,143,36,241,17,185,171,58,34,86,169,17,185,85,61,50,162,2,43,126,66,36,67,67,99,29,33,139,184,204,11,13,57,113,132,12,209,51,61,52,164,198,
		17,49,187,10,145,121,166,163,69,198,146,181,141,132,245,145,192,254,200,88,210,181,145,176,63,84,244,106,173,161,161,181,13,17,89,166,133,200,46,243,34,66,100,158,233,136,33,161,171,208,
		144,195,101,164,186,12,13,181,117,104,168,149,66,69,181,42,21,26,170,212,208,88,91,41,161,161,51,143,200,169,169,145,184,127,68,174,178,130,144,104,1,
	};
	const unsigned char sf6[]=//small fonts size 6
	{
		145,168,27,17,249,8,149,189,189,69,198,244,249,202,16,37,73,70,168,164,180,52,18,114,68,204,42,35,98,169,22,17,249,204,140,174,88,200,33,137,99,17,35,82,122,69,166,212,86,17,
		41,87,35,83,163,242,200,212,40,46,50,69,247,140,76,61,231,34,83,202,170,200,212,41,41,50,165,170,138,76,169,167,66,49,13,5,117,72,198,76,73,28,146,201,34,82,25,69,204,24,
		169,237,5,38,66,101,150,159,17,170,203,203,139,80,89,70,90,132,234,50,243,34,83,207,242,200,212,179,36,66,229,209,233,17,170,204,207,140,68,253,200,20,201,46,66,149,53,149,145,169,
		36,121,164,42,238,26,99,132,42,183,51,35,84,150,153,22,161,186,188,136,8,149,101,150,70,168,46,47,51,66,229,97,120,145,169,75,82,132,42,51,211,34,83,109,165,72,85,172,85,169,
		200,84,171,54,50,213,74,138,76,157,202,35,98,87,29,145,170,53,34,182,234,145,17,21,71,241,17,145,12,205,184,142,76,37,235,66,51,142,35,83,164,235,208,140,186,136,148,93,161,33,
		215,140,76,37,107,35,81,61,18,214,71,166,146,174,145,168,31,170,121,181,134,102,214,134,102,222,135,134,214,37,52,228,154,33,153,11,52,99,45,35,212,25,154,105,29,154,105,21,162,201,
		252,208,76,213,208,80,165,132,102,174,71,196,180,52,18,245,35,98,165,21,25,241,0,
	};
}
#endif
void			uncompress_bitmap_v4(const unsigned char *data, int data_size, char fontH, const char *widths, int *&rgb2, int &w2, int &h2)
{
	const int red=0xFF0000FF, green=0xFF00FF00, blue=0xFFFF0000;//swapped in WinAPI bitmaps
	const int size=txw*txh;
	w2=txw, h2=txh;
	rgb2=(int*)malloc(size<<2);
	for(int k=0;k<size;++k)
		rgb2[k]=blue;
	int nbits=data_size<<3;
	std::vector<bool> cdata(nbits);
	for(int k=0;k<nbits;k+=8)
	{
		for(int bit=0;bit<8;++bit)
			cdata[k+bit]=data[k>>3]>>bit&1;
	}
	{
		int c=' ';
		int xstart=(c&7)<<4, ystart=(c>>3)<<4;
		int bkwidth=widths[c];
		for(int ky=0;ky<fontH;++ky)
		{
			for(int kx=0;kx<bkwidth;++kx)
			{
				rgb2[(ystart+ky)<<7|(xstart+kx)]=green;
			}
		}
	}
	for(int c='!', pos=0;c<127;++c)//95 printable characters, space is all-bk
	{
		int xstart=(c&7)<<4, ystart=(c>>3)<<4;
		int bkwidth=widths[c];
		for(int ky=0;ky<fontH;++ky)
		{
			for(int kx=0;kx<bkwidth;++kx)
			{
				rgb2[(ystart+ky)<<7|(xstart+kx)]=green;
			}
		}
		int	xoffset=(int)cdata[pos+2]<<2|(int)cdata[pos+1]<<1|(int)cdata[pos],
			yoffset=(int)cdata[pos+6]<<3|(int)cdata[pos+5]<<2|(int)cdata[pos+4]<<1|(int)cdata[pos+3],
			width=(int)cdata[pos+10]<<3|(int)cdata[pos+9]<<2|(int)cdata[pos+8]<<1|(int)cdata[pos+7],
			height=(int)cdata[pos+14]<<3|(int)cdata[pos+13]<<2|(int)cdata[pos+12]<<1|(int)cdata[pos+11];
		xstart+=xoffset, ystart+=yoffset;
		pos+=15;
		for(int ky=0;ky<height;++ky)
		{
			for(int kx=0;kx<width;++kx)
			{
				if(cdata[pos])
					rgb2[(ystart+ky)<<7|(xstart+kx)]=red;
				++pos;
			}
		}
	}
}
void			LoadTexture(const unsigned char *data, int data_size, char fontH, const char *widths, unsigned &tx_id)
{
	int w2, h2;
	int *rgb2=nullptr;
	uncompress_bitmap_v4(data, data_size, fontH, widths, rgb2, w2, h2);
	LoadTexture(rgb2, w2, h2, tx_id, false);
	free(rgb2);
}

namespace		GL2_Text//font & textures
{
	unsigned		program=0;
	int				a_coord2d=-1, u_mytexture=-1, u_txtColor=-1, u_bkColor=-1, u_isTexture=-1;
	unsigned		buffer=0;
}

unsigned		tx_id_sf10=0, tx_id_sf8=0, tx_id_sf7=0, tx_id_sf6=0,
				font_tx_id=0;
float			sf10_txcoords[128<<2]={0}, sf8_txcoords[128<<2]={0}, sf7_txcoords[128<<2]={0}, sf6_txcoords[128<<2]={0},//txx1 txx2 txy1 txy2
				*text_txcoords=nullptr;
const char		*font_widths=nullptr;
int				gl_tabWidth=0;
int				gl_font_size=0,
				pixel_x=2, pixel_y=2, gl_fontH=sf10_height,
				txtColor=0xCF0000FF,//0xFF000000 0x7F000000 0x7F00FF00	//WinAPI DIB: 0xAARRGGBB, WinAPI functions: 0xAABBGGRR, OpenGL: 0xAABBGGRR
				bkColor=0x3FFF0000;//0x00FFFFFF 0x707F7F7F 0x3F23B3BA
void			gl_setTextColor(int color)
{
	gl_setProgram(GL2_Text::program);
	send_color(GL2_Text::u_txtColor, color);
}
void			gl_setBkColor(int color)
{
	gl_setProgram(GL2_Text::program);
	send_color(GL2_Text::u_bkColor, color);
}
int				gl_setTextSize(int size)//{1, ..., 10}
{
	switch(size=clamp(1, size, 10))
	{
	case  1:font_tx_id=tx_id_sf6,	text_txcoords=sf6_txcoords,		font_widths=sf6_widths,		pixel_x=pixel_y=1,		gl_fontH= 8, gl_tabWidth=32;	break;
	case  2:font_tx_id=tx_id_sf7,	text_txcoords=sf7_txcoords,		font_widths=sf7_widths,		pixel_x=pixel_y=1,		gl_fontH=10, gl_tabWidth=40;	break;
	case  3:font_tx_id=tx_id_sf8,	text_txcoords=sf8_txcoords,		font_widths=sf8_widths,		pixel_x=pixel_y=1,		gl_fontH=12, gl_tabWidth=48;	break;
	case  4:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=1,		gl_fontH=16, gl_tabWidth=64;	break;
	case  5:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=2,		gl_fontH=16, gl_tabWidth=136;	break;
	case  6:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=3,		gl_fontH=16, gl_tabWidth=200;	break;
	case  7:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=pixel_y=4,		gl_fontH=16, gl_tabWidth=264;	break;
	case  8:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=5, pixel_y=6,	gl_fontH=16, gl_tabWidth=328;	break;
	case  9:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=5, pixel_y=7,	gl_fontH=16, gl_tabWidth=328;	break;
	case 10:font_tx_id=tx_id_sf10,	text_txcoords=sf10_txcoords,	font_widths=sf10_widths,	pixel_x=5, pixel_y=8,	gl_fontH=16, gl_tabWidth=328;	break;
	}
	return size;
}
inline int		gl_getFontHeight(){return pixel_y*gl_fontH;}
inline int		gl_getCharWidth(char c){return pixel_x*font_widths[c];}
std::vector<float> vrtx;
int				print_array(int x, int y, const char *msg, int msg_length, int tab_origin)//stutters on RX 580?, tab origin		6.4fps
{
	gl_setProgram(GL2_Text::program);				check(__LINE__);
	glUniform2f(GL2_Text::u_isTexture, 0, 1);		check(__LINE__);
	int msg_width=0;
	float _2_w=2.f/w, _2_h=2.f/h;
	select_texture(font_tx_id, GL2_Text::u_mytexture);
	vrtx.resize(msg_length<<4);//vx, vy, txx, txy		6.5fps
	float rect[4], *txc;
	int width, idx;
	int fontH_px=gl_fontH*pixel_y;
	for(int k=0;k<msg_length;++k)
	{
		char c=msg[k];
		width=0;
		if(c=='\t')
			width=gl_tabWidth-(msg_width-tab_origin)%gl_tabWidth, c=' ';
		else if(c>=32&&c<0xFF&&(width=font_widths[c]))
			width*=pixel_x;
		rect[0]=(x+msg_width		)*_2_w-1, rect[1]=1- y			*_2_h;
		rect[2]=(x+msg_width+width	)*_2_w-1, rect[3]=1-(y+fontH_px)*_2_h;//y2<y1
		idx=k<<4;
		txc=text_txcoords+(c<<2);
		vrtx[idx   ]=rect[0], vrtx[idx+ 1]=rect[1],		vrtx[idx+ 2]=txc[0], vrtx[idx+ 3]=txc[2];//top left
		vrtx[idx+ 4]=rect[0], vrtx[idx+ 5]=rect[3],		vrtx[idx+ 6]=txc[0], vrtx[idx+ 7]=txc[3];//bottom left
		vrtx[idx+ 8]=rect[2], vrtx[idx+ 9]=rect[3],		vrtx[idx+10]=txc[1], vrtx[idx+11]=txc[3];//bottom right
		vrtx[idx+12]=rect[2], vrtx[idx+13]=rect[1],		vrtx[idx+14]=txc[1], vrtx[idx+15]=txc[2];//top right

		msg_width+=width;
	}
	glBindBuffer(GL_ARRAY_BUFFER, GL2_Text::buffer);									check(__LINE__);
	glBufferData(GL_ARRAY_BUFFER, msg_length<<6, &vrtx[0], GL_STATIC_DRAW);				check(__LINE__);//set vertices & texcoords
	glVertexAttribPointer(GL2_Text::a_coord2d, 4, GL_FLOAT, GL_FALSE, 4<<2, (void*)0);	check(__LINE__);
	glEnableVertexAttribArray(GL2_Text::a_coord2d);		check(__LINE__);
	glDrawArrays(GL_QUADS, 0, msg_length<<2);			check(__LINE__);//draw the quads: 4 vertices per character quad
	glDisableVertexAttribArray(GL2_Text::a_coord2d);	check(__LINE__);
	return msg_width;
}
int				gl_print(int x, int y, int tab_origin, const char *format, ...)
{
	int msg_length=vsnprintf_s(g_buf, g_buf_size, format, (char*)(&format+1));
	return print_array(x, y, g_buf, msg_length, tab_origin);//38fps 41fps	48fps		undebugged 6.4fps
//	return print_array(x, y, g_buf, msg_length);//7.6fps
}
void			display_texture(int x1, int x2, int y1, int y2, int *rgb, int txw, int txh, unsigned char alpha=0xFF)
{
	static unsigned tx_id=0;
	float _2_w=2.f/w, _2_h=2.f/h;
	float rect[]=
	{
		x1*_2_w-1, 1-y1*_2_h,
		x2*_2_w-1, 1-y2*_2_h//y2<y1
	};
	float vrtx[]=
	{
		rect[0], rect[1],		0, 0,//top left
		rect[0], rect[3],		0, 1,//bottom left
		rect[2], rect[3],		1, 1,//bottom right
		rect[2], rect[1],		1, 0 //top right
	};
	if(rgb)
	{
	#define	NPOT_ATIX2300_FIX
		int *rgb2, w2, h2;
#ifdef NPOT_ATIX2300_FIX
		bool expand;
		int logw=floor_log2(txw), logh=floor_log2(txh);
		if(expand=glMajorVer<3&&(txw>1<<logw||txh>1<<logh))
		{
			w2=txw>1<<logw?1<<(logw+1):txw;
			h2=txh>1<<logh?1<<(logh+1):txh;
			int size=w2*h2;
			rgb2=(int*)malloc(size<<2);
			memset(rgb2, 0, size<<2);
			for(int ky=0;ky<txh;++ky)
				memcpy(rgb2+w2*ky, rgb+txw*ky, txw<<2);
		//	memcpy(rgb2, rgb, size<<2);
			float nw=(float)txw/w2, nh=(float)txh/h2;
			vrtx[ 2]=0,		vrtx[ 3]=0;
			vrtx[ 6]=0,		vrtx[ 7]=nh;
			vrtx[10]=nw,	vrtx[11]=nh;
			vrtx[14]=nw,	vrtx[15]=0;
		}
		else
#endif
			rgb2=rgb, w2=txw, h2=txh;
		gl_setProgram(GL2_Text::program);				check(__LINE__);//select Text program
		glUniform2f(GL2_Text::u_isTexture, 1, 0);		check(__LINE__);//send isTexture
	//	glUniform1i(GL2_Text::u_isTexture, true);		check(__LINE__);
		send_color(GL2_Text::u_bkColor, alpha<<24);		check(__LINE__);//send apha
		if(!tx_id)
			{glGenTextures(1, &tx_id);														check(__LINE__);}//generate texture id once
		glBindTexture(GL_TEXTURE_2D, tx_id);												check(__LINE__);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);					check(__LINE__);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);					check(__LINE__);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w2, h2, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb2);check(__LINE__);//send bitmap to GPU

		select_texture(tx_id, GL2_Text::u_mytexture);
		glBindBuffer(GL_ARRAY_BUFFER, GL2_Text::buffer);									check(__LINE__);
		glBufferData(GL_ARRAY_BUFFER, 16<<2, vrtx, GL_STATIC_DRAW);							check(__LINE__);//send vertices & texcoords

		glVertexAttribPointer(GL2_Text::a_coord2d, 4, GL_FLOAT, GL_FALSE, 4<<2, (void*)0);	check(__LINE__);//select vertices & texcoord

		glEnableVertexAttribArray(GL2_Text::a_coord2d);		check(__LINE__);
		glDrawArrays(GL_QUADS, 0, 4);						check(__LINE__);//draw the quad
		glDisableVertexAttribArray(GL2_Text::a_coord2d);	check(__LINE__);
#ifdef NPOT_ATIX2300_FIX
		if(expand)
			free(rgb2);
#endif
	}
}

int			getBkMode(){return 1+(((unsigned char*)&bkColor)[3]==0xFF);}
int			setBkMode(int mode)
{
	gl_setProgram(GL2_Text::program);
	unsigned char alpha=((unsigned char*)&bkColor)[3];
	((unsigned char*)&bkColor)[3]=-(mode==2);
	send_color(GL2_Text::u_bkColor, bkColor);
	return 1+(alpha==0xFF);
}
int			getBkColor(){return bkColor&0xFFFFFF;}
int			setBkColor(int color)
{
	gl_setProgram(GL2_Text::program);
	int oldColor=bkColor&0xFFFFFF;
	((unsigned char*)&bkColor)[0]=((unsigned char*)&color)[0];
	((unsigned char*)&bkColor)[1]=((unsigned char*)&color)[1];
	((unsigned char*)&bkColor)[2]=((unsigned char*)&color)[2];
	send_color(GL2_Text::u_bkColor, bkColor);
	return oldColor;
}
int			getTextColor(){return txtColor&0xFFFFFF;}
int			setTextColor(int color)
{
	gl_setProgram(GL2_Text::program);
	int oldColor=txtColor;
	txtColor=0xFF000000|0xFFFFFF&color;
	send_color(GL2_Text::u_txtColor, txtColor);
	return oldColor;
}
struct Text
{
	int x, y;
	int bkMode;
	bool enable_color;	int color;
	std::string str;
	Text(int x, int y, int bkMode, std::string &str)			:x(x), y(y), bkMode(bkMode), enable_color(false), str(str){}
	Text(int x, int y, int bkMode, const char *a)				:x(x), y(y), bkMode(bkMode), enable_color(false), str(a){}
	Text(int x, int y, int bkMode, int color, std::string &str)	:x(x), y(y), bkMode(bkMode), enable_color(true), color(color), str(str){}
	Text(int x, int y, int bkMode, int color, const char *a)	:x(x), y(y), bkMode(bkMode), enable_color(true), color(color), str(a){}
};
std::map<float, std::list<Text>> LOL_text;//closest to farthest
void label(float x, float y, float z, char const *format, ...)
{
	vec3 cp;
	cam.world2cam(vec3(x, y, z), cp);
	if(cp.z>0)
	{
		vec2 s;
		cam.cam2screen(cp, s);
		if(abs(s.x)<1e6&&abs(s.y)<1e6&&vsprintf_s(g_buf, g_buf_size, format, (char*)(&format+1))>0)
			LOL_text[1/cp.z].push_back(Text(int(s.x)-(s.x<0), int(s.y)-(s.y<0), TRANSPARENT, g_buf));//-0.5 truncated as 0
	}
}
void text_show()//show the text that was printed in 3d space
{
	for(auto it=LOL_text.begin();it!=LOL_text.end();++it)
	{
		auto &LOL=*it;
		for(auto it2=LOL.second.begin();it2!=LOL.second.end();++it2)
		{
			auto &LOL_2=*it2;
			//setBkMode(LOL_2.bkMode);
			//if(LOL_2.enable_color)
			//	setTextColor(LOL_2.color);
			gl_print(LOL_2.x, LOL_2.y, 0, LOL_2.str.c_str(), LOL_2.str.size());
			//if(LOL_2.enable_color)
			//	setTextColor(0);
		}
	}
	LOL_text.clear();
}

void			push_vec(std::vector<float> &vertices, vec3 const &v)
{
	vertices.reserve(vertices.size()+3);
	vertices.push_back(v.x), vertices.push_back(v.y), vertices.push_back(v.z);
}
void			push_vec(std::vector<float> &vertices, vec2 const &v)
{
	vertices.reserve(vertices.size()+2);
	vertices.push_back(v.x), vertices.push_back(v.y);
}
void			generate_oh_face(std::vector<float> &vertices, std::vector<int> &indices, vec3 const &c, vec3 const &p1, vec3 const &p2, vec3 const &p3, unsigned resolution, bool normals, bool texcoords)
{
	float inv_res=1.f/resolution;
	int stride=3+3*normals+2*texcoords;
	vec3 u12=(p2-p1)*inv_res, u23=(p3-p2)*inv_res;
	vec3 n;
	int idx=vertices.size()/stride;
	push_vec(vertices, p1);//push top point
	if(normals)
	{
		n=p1-c;
		push_vec(vertices, n/n.magnitude());
	}
	if(texcoords)
		push_vec(vertices, vec2(0, 0));
	for(int k=1;k<=(int)resolution;++k)//for each row
	{
		vec3 start=p1+(float)k*u12;
		push_vec(vertices, start);//push row start point
		if(normals)
		{
			n=start-c;
			push_vec(vertices, n/n.magnitude());
		}
		if(texcoords)
			push_vec(vertices, vec2(0, 0));
		for(int k2=1;k2<=k;++k2, ++idx)//for each point
		{
			vec3 np=start+(float)k2*u23;
			push_vec(vertices, np);
			if(normals)
			{
				n=np-c;
				push_vec(vertices, n/n.magnitude());
			}
			if(texcoords)
				push_vec(vertices, vec2(0, 0));
			indices.push_back(idx);
			indices.push_back(vertices.size()/stride-2);
			indices.push_back(vertices.size()/stride-1);
			if(k2<k)
			{
				indices.push_back(vertices.size()/stride-1);
				indices.push_back(idx+1);
				indices.push_back(idx);
			}
		}
	}
}
void			normalize_about(std::vector<float> &vertices, vec3 const &center, float radius, bool normals, bool texcoords)
{
	int stride=3+3*normals+2*texcoords;
	for(int k=0, kEnd=vertices.size();k<kEnd;k+=stride)
	{
		vec3 p(&vertices[k]);
		auto delta=p-center;
		p=center+delta*(radius/delta.magnitude());
		memcpy(&vertices[k], &p.x, 3<<2);
	}
}
void			generate_sphere(std::vector<float> &vertices, std::vector<int> &indices, vec3 const &center, float radius, unsigned resolution, bool normals, bool texcoords)
{
	vec3
		up(0, 0, radius), north(0, radius, 0), east(radius, 0, 0),
		vu=center+up, vn=center+north, ve=center+east,
		vd=center-up, vs=center-north, vw=center-east;
	generate_oh_face(vertices, indices, center, vu, vn, ve, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vu, vs, ve, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vu, vn, vw, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vu, vs, vw, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vd, vn, ve, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vd, vs, ve, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vd, vn, vw, resolution, normals, texcoords);
	generate_oh_face(vertices, indices, center, vd, vs, vw, resolution, normals, texcoords);

	//vbo_to_clipboard(vertices, indices, true, true);
	normalize_about(vertices, center, radius, normals, texcoords);
	//vbo_to_clipboard(vertices, indices, true, true);
}
int				s_vcount=0;

unsigned		vbo_cube_vertices, vbo_cube_colors, ibo_cube_elements, vbo_cube_texcoords;
float			cube_vertices[]=
{
	-1, -1,  1,//front
	 1, -1,  1,
	 1,  1,  1,
	-1,  1,  1,
	
	-1, -1, -1,//back
	 1, -1, -1,
	 1,  1, -1,
	-1,  1, -1
};
float			cube_colors[]=
{
	1, 0, 0,//front colors
	0, 1, 0,
	0, 0, 1,
	1, 1, 1,
    
	1, 0, 0,//back colors
	0, 1, 0,
	0, 0, 1,
	1, 1, 1
};
unsigned short	cube_elements[]=
{
	0, 1, 2,	2, 3, 0,//front
	1, 5, 6,	6, 2, 1,//right
	7, 6, 5,	5, 4, 7,//back
	4, 0, 3,	3, 7, 4,//left
	4, 5, 1,	1, 0, 4,//bottom
	3, 2, 6,	6, 7, 3 //top
};
float			cube_texcoords[2*4*6]=
{
	0, 1,//front
	1, 1,
	1, 0,
	0, 0,
};
std::wstring	wider(const char *a)
{
	auto len=strlen(a);
	auto buf=new wchar_t[len+1];
	swprintf_s(buf, len+1, L"%S", a);
	std::wstring str=buf;
	delete[] buf;
	return str;
}
void			openImage(const char *filename, int *&rgb, int &width, int &height, bool swap_rb)
{
	int size=0;
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR hgdiplusToken;
	Gdiplus::GdiplusStartup(&hgdiplusToken, &gdiplusStartupInput, 0);
	{
		Gdiplus::Bitmap bitmap(wider(filename).c_str());
		height=bitmap.GetHeight(), width=bitmap.GetWidth();
		size=height*width;
		rgb=(int*)malloc(size<<2);
		Gdiplus::Rect rt(0, 0, width, height);
		Gdiplus::BitmapData data;
		bitmap.LockBits(&rt, Gdiplus::ImageLockModeRead, PixelFormat32bppARGB, &data);
		memcpy(rgb, data.Scan0, size<<2);
		bitmap.UnlockBits(&data);
	}
	if(swap_rb)
		for(int k=0;k<size;++k)//swap red and blue channels
		{
			auto p=(unsigned char*)(rgb+k);
			unsigned char temp=p[0]; p[0]=p[2], p[2]=temp;
		}
	Gdiplus::GdiplusShutdown(hgdiplusToken);
}
unsigned		LoadTexture(const char *filename)
{
    unsigned texture;
	
	int *rgb;
	int height, width;
	openImage(filename, rgb, width, height, true);

	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgb);
	check(__LINE__);
	
	free(rgb);
    return texture;
}
unsigned		make_gpu_buffer(unsigned target, const void *pointer, int size_bytes)
{
	unsigned buffer_id=0;
	glGenBuffers(1, &buffer_id);
	glBindBuffer(target, buffer_id);
	glBufferData(target, size_bytes, pointer, GL_STATIC_DRAW);
	return buffer_id;
}
void			calculate_text_txcoords(const char *font_widths, int gl_fontH, float *txcoord)
{
	for(int c=0, idx=0;c<128;++c, idx+=4)
	{
		int width=font_widths[c];
		int xpos=(c&0x7)<<4, ypos=(c>>3&0xF)<<4;
		txcoord[idx  ]=xpos*inv128, txcoord[idx+1]=(xpos+width)*inv128;
		txcoord[idx+2]=ypos*inv256, txcoord[idx+3]=(ypos+gl_fontH)*inv256;
	}
}
struct			LinearGL:private Engine
{
	HGLRC__		*hRC;
	int			initiate()
	{
		prof_start();//

		//HDC dc=GetDC(0);
		//DescribePixelFormat(dc, 0, 0, 0);
		//ReleaseDC(NULL, dc);

		tagPIXELFORMATDESCRIPTOR pfd={sizeof(tagPIXELFORMATDESCRIPTOR), 1, PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER, PFD_TYPE_RGBA, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, PFD_MAIN_PLANE, 0, 0, 0, 0};
		int PixelFormat=ChoosePixelFormat(ghDC, &pfd);
		prof_add("ChoosePixelFormat");//either this takes 60ms
		SetPixelFormat(ghDC, PixelFormat, &pfd);
		prof_add("SetPixelFormat");//...or this 60ms
		hRC=wglCreateContext(ghDC);
		prof_add("wglCreateContext");//37.5ms
		wglMakeCurrent(ghDC, hRC);
		prof_add("wglMakeCurrent");
		if(h<=0)
			h=1;
		changeFov();
		prof_add("glViewport");
		glShadeModel(GL_SMOOTH);
		prof_add("glShadeModel");
		glClearColor(1, 1, 1, 1);
	//	glClearColor(0, 0, 0, 1);
		prof_add("glClearColor");
		glClearDepth(1);
		prof_add("glClearDepth");
		glEnable(GL_DEPTH_TEST);
		prof_add("glEnable");
		glDepthFunc(GL_LEQUAL);
		prof_add("glDepthFunc");
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		prof_add("glHint");

		glGetIntegerv(GL_MAJOR_VERSION, &glMajorVer);
		glGetIntegerv(GL_MINOR_VERSION, &glMinorVer);
		GLversion=glGetString(GL_VERSION);
		if(API_not_loaded)
		{
		//	glGenVertexArrays=(decltype(glGenVertexArrays))wglGetProcAddress("glGenVertexArrays");//OpenGL 3.0
		//	glDeleteVertexArrays=(decltype(glDeleteVertexArrays))wglGetProcAddress("glDeleteVertexArrays");//OpenGL 3.0
			glBindVertexArray=(decltype(glBindVertexArray))wglGetProcAddress("glBindVertexArray");
			glGenBuffers=(decltype(glGenBuffers))wglGetProcAddress("glGenBuffers");
			glBindBuffer=(decltype(glBindBuffer))wglGetProcAddress("glBindBuffer");
			glBufferData=(decltype(glBufferData))wglGetProcAddress("glBufferData");
			glBufferSubData=(decltype(glBufferSubData))wglGetProcAddress("glBufferSubData");
			glEnableVertexAttribArray=(decltype(glEnableVertexAttribArray))wglGetProcAddress("glEnableVertexAttribArray");
			glVertexAttribPointer=(decltype(glVertexAttribPointer))wglGetProcAddress("glVertexAttribPointer");
			glDisableVertexAttribArray=(decltype(glDisableVertexAttribArray))wglGetProcAddress("glDisableVertexAttribArray");
			glCreateShader=(decltype(glCreateShader))wglGetProcAddress("glCreateShader");
			glShaderSource=(decltype(glShaderSource))wglGetProcAddress("glShaderSource");
			glCompileShader=(decltype(glCompileShader))wglGetProcAddress("glCompileShader");
			glGetShaderiv=(decltype(glGetShaderiv))wglGetProcAddress("glGetShaderiv");
			glGetShaderInfoLog=(decltype(glGetShaderInfoLog))wglGetProcAddress("glGetShaderInfoLog");
			glCreateProgram=(decltype(glCreateProgram))wglGetProcAddress("glCreateProgram");
			glAttachShader=(decltype(glAttachShader))wglGetProcAddress("glAttachShader");
			glLinkProgram=(decltype(glLinkProgram))wglGetProcAddress("glLinkProgram");
			glGetProgramiv=(decltype(glGetProgramiv))wglGetProcAddress("glGetProgramiv");
			glGetProgramInfoLog=(decltype(glGetProgramInfoLog))wglGetProcAddress("glGetProgramInfoLog");
			glDetachShader=(decltype(glDetachShader))wglGetProcAddress("glDetachShader");
			glDeleteShader=(decltype(glDeleteShader))wglGetProcAddress("glDeleteShader");
			glUseProgram=(decltype(glUseProgram))wglGetProcAddress("glUseProgram");
			glGetAttribLocation=(decltype(glGetAttribLocation))wglGetProcAddress("glGetAttribLocation");
			glDeleteProgram=(decltype(glDeleteProgram))wglGetProcAddress("glDeleteProgram");
			glDeleteBuffers=(decltype(glDeleteBuffers))wglGetProcAddress("glDeleteBuffers");
			glGetUniformLocation=(decltype(glGetUniformLocation))wglGetProcAddress("glGetUniformLocation");
			glUniform1f=(decltype(glUniform1f))wglGetProcAddress("glUniform1f");
			glUniformMatrix3fv=(decltype(glUniformMatrix3fv))wglGetProcAddress("glUniformMatrix3fv");
			glUniformMatrix4fv=(decltype(glUniformMatrix4fv))wglGetProcAddress("glUniformMatrix4fv");
			glGetBufferParameteriv=(decltype(glGetBufferParameteriv))wglGetProcAddress("glGetBufferParameteriv");
			glActiveTexture=(decltype(glActiveTexture))wglGetProcAddress("glActiveTexture");
			glUniform1i=(decltype(glUniform1i))wglGetProcAddress("glUniform1i");
			glUniform2f=(decltype(glUniform2f))wglGetProcAddress("glUniform2f");
			glUniform3f=(decltype(glUniform3f))wglGetProcAddress("glUniform3f");
			glUniform3fv=(decltype(glUniform3fv))wglGetProcAddress("glUniform3fv");
			glUniform4f=(decltype(glUniform4f))wglGetProcAddress("glUniform4f");
			API_not_loaded=false;
		}
		prof_add("Load API");
#if 1
		vbo_cube_vertices	=make_gpu_buffer(GL_ELEMENT_ARRAY_BUFFER, cube_vertices, sizeof(cube_vertices));	check(__LINE__);
		vbo_cube_colors		=make_gpu_buffer(GL_ELEMENT_ARRAY_BUFFER, cube_colors, sizeof(cube_colors));		check(__LINE__);
		ibo_cube_elements	=make_gpu_buffer(GL_ELEMENT_ARRAY_BUFFER, cube_elements, sizeof(cube_elements));	check(__LINE__);
		for(int i=1;i<6;++i)
			memcpy(&cube_texcoords[i*4*2], &cube_texcoords[0], 2*4*sizeof(float));
		vbo_cube_texcoords	=make_gpu_buffer(GL_ARRAY_BUFFER, cube_texcoords, sizeof(cube_texcoords));			check(__LINE__);

		glGenBuffers(1, &GL2_2D::vertex_buffer);
	//	GL2_2D::vertex_buffer=make_gpu_buffer(GL_ARRAY_BUFFER, 0, 2<<2);//dummy size

		//{//generate & send sphere
		//	std::vector<float> vertices;
		//	std::vector<int> indices;
		//	generate_sphere(vertices, indices, vec3(), 5, 9, true, false);
		////	generate_sphere(vertices, indices, vec3(), 5, 3, true, false);//10*8=80 vertices
		////	generate_sphere(vertices, indices, vec3(), 5, 1, true, false);//3*8=24 vertices octahedron
		//	s_vcount=indices.size();

		//	glGenBuffers(1, &sphereVBO);
		//	glGenBuffers(1, &sphereEBO);
		//	glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);												check(__LINE__);
		//	glBufferData(GL_ARRAY_BUFFER, vertices.size()<<2, &vertices[0], GL_STATIC_DRAW);		check(__LINE__);
		//	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereEBO);										check(__LINE__);
		//	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()<<2, &indices[0], GL_STATIC_DRAW);	check(__LINE__);
		//}

		//glGenBuffers(1, &GL2_BB3D::buffer);

		glGenBuffers(1, &GL2_3D::vertex_buffer);			check(__LINE__);
		//float test_triangle[]=
		//{
		//	-0.5f, -0.5f, 0,
		//	0.5f, -0.5f, 0,
		//	0, 0.5f, 0
		//};
	//	GL2_3D::vertex_buffer=make_gpu_buffer(GL_ARRAY_BUFFER, 0, 3<<2);//dummy sizes
		float texcoords[]=
		{
			0, 0,
			0, 1,
			1, 1,
			1, 0
		};
		GL2_3D::texcoord_buffer=make_gpu_buffer(GL_ARRAY_BUFFER, texcoords, sizeof(texcoords));
		
		calculate_text_txcoords(sf10_widths, sf10_height, sf10_txcoords);
		calculate_text_txcoords(sf8_widths, sf8_height, sf8_txcoords);
		calculate_text_txcoords(sf7_widths, sf7_height, sf7_txcoords);
		calculate_text_txcoords(sf6_widths, sf6_height, sf6_txcoords);
		GL2_Text::buffer=make_gpu_buffer(GL_ARRAY_BUFFER, 0, 8<<2);		check(__LINE__);//dummy size
		prof_add("Alloc bufs");

		ShaderVar ex_attr[]=					//Learning OpenGL Example
		{
			{&attribute_coord3d, "coord3d", __LINE__},
			{&attribute_texcoord, "texcoord", __LINE__}
		};
		ShaderVar ex_unif[]=
		{
			{&uniform_mvp, "mvp", __LINE__},
			{&uniform_texture, "mytexture", __LINE__},
			{&uniform_fade, "fade", __LINE__}
		};
		programID=LoadShaders(
			"#version 120\n"
			"attribute vec3 coord3d;\n"
			"attribute vec2 texcoord;\n"
			"varying vec2 f_texcoord;\n"
			"uniform mat4 mvp;\n"
			"void main()\n"
			"{\n"
			"    gl_Position = mvp*vec4(coord3d, 1.0);\n"
			"    f_texcoord = texcoord;\n"
			"}",
			"#version 120\n"
			"varying vec2 f_texcoord;\n"
			"uniform sampler2D mytexture;\n"
			"uniform float fade;\n"
			"void main()\n"
			"{\n"
			"	 vec4 color=texture2D(mytexture, f_texcoord);\n"
			"    gl_FragColor=vec4(color.rgb, fade);\n"
			"}",
			ex_attr, sizeof(ex_attr)/sizeof(ShaderVar), ex_unif, sizeof(ex_unif)/sizeof(ShaderVar));
		if(!programID)
			error(__LINE__);
		
		ShaderVar _2d_attr[]=				//2D program
		{
			{&GL2_2D::a_vertices, "a_vertex", __LINE__}
		};
		ShaderVar _2d_unif[]=
		{
			{&GL2_2D::u_color, "u_color", __LINE__}
		};
		GL2_2D::program=LoadShaders(
			"#version 120\n"
			"attribute vec2 a_vertex;\n"
			"void main()\n"
			"{\n"
			"    gl_Position=vec4(a_vertex, 0., 1.);\n"
			"}",
			"#version 120\n"
			"uniform vec4 u_color;\n"
			"void main()\n"
			"{\n"
			"    gl_FragColor=u_color;\n"
			"}",
			_2d_attr, sizeof(_2d_attr)/sizeof(ShaderVar), _2d_unif, sizeof(_2d_unif)/sizeof(ShaderVar));
		if(!GL2_2D::program)
			error(__LINE__);
		
		ShaderVar l3d_attr[]=						//Light 3D program
		{
			{&GL2_L3D::a_vertices, "a_vertex", __LINE__},
			{&GL2_L3D::a_normals, "a_normal", __LINE__},
		//	{&GL2_L3D::a_texcoord, "a_texcoord", __LINE__}
		};
		ShaderVar l3d_unif[]=
		{
			{&GL2_L3D::u_vpmatrix, "u_vpmatrix", __LINE__},
			{&GL2_L3D::u_modelmatrix, "u_modelmatrix", __LINE__},
			{&GL2_L3D::u_normalmatrix, "u_normalmatrix", __LINE__},
		//	{&GL2_L3D::u_pointSize, "u_pointSize", __LINE__},
		//	{&GL2_L3D::u_texture, "u_texture", __LINE__},
		//	{&GL2_L3D::u_isTextured, "u_isTextured", __LINE__},
			{&GL2_L3D::u_objectcolor, "u_objectcolor", __LINE__},
			{&GL2_L3D::u_lightcolor, "u_lightcolor", __LINE__},
			{&GL2_L3D::u_lightpos, "u_lightpos", __LINE__},
			{&GL2_L3D::u_viewpos, "u_viewpos", __LINE__},
		};
		GL2_L3D::program=LoadShaders(
			"#version 120\n"
			"uniform mat4 u_vpmatrix, u_modelmatrix;\n"
			"uniform mat3 u_normalmatrix;\n"
			"attribute vec3 a_vertex;\n"
			"attribute vec3 a_normal;\n"
		//	"attribute vec2 a_texcoord;\n"
			"varying vec3 v_fragpos;\n"
			"varying vec3 v_normal;\n"
		//	"varying vec2 v_texcoord;\n"
		//	"uniform float u_pointSize;\n"//can't get u_pointSize, -1
			"varying vec4 v_glposition;\n"
			"void main()\n"
			"{\n"
			"    vec4 fullpos=vec4(a_vertex, 1.);\n"
			"    gl_Position=u_vpmatrix*fullpos;\n"
			"    v_glposition=gl_Position;\n"
			"    gl_Position.z=0.;\n"
			"    v_fragpos=vec3(u_modelmatrix*fullpos);\n"
			"    v_normal=u_normalmatrix*a_normal;\n"
		//	"    v_texcoord=a_texcoord;\n"
		//	"    gl_PointSize=u_pointSize;\n"
			"    gl_PointSize=20.;\n"
			"}",
			"#version 120\n"
			"varying vec3 v_fragpos;\n"
			"varying vec3 v_normal;\n"
			"varying vec4 v_glposition;\n"
		//	"varying vec2 v_texcoord;\n"
		//	"uniform bool u_isTextured;\n"
			"uniform vec4 u_objectcolor;\n"
			"uniform vec3 u_lightcolor, u_lightpos, u_viewpos;\n"
		//	"uniform sampler2D u_texture;\n"
			"void main()\n"
			"{\n"
			"    vec3 normal=normalize(v_normal);\n"
			"    vec3 lightdir=normalize(u_lightpos-v_fragpos);\n"
				
			"    float specularstrength=0.5;\n"
			"    vec3 viewdir=normalize(u_viewpos-v_fragpos), reflectdir=reflect(-lightdir, normal);\n"
			"    vec3 specular=specularstrength*u_lightcolor*pow(max(dot(viewdir, reflectdir), 0.), 32);\n"

			"    vec3 diffuse=max(dot(normal, lightdir), 0.)*u_lightcolor;\n"
			"    gl_FragColor=vec4((0.1*u_lightcolor+diffuse+specular)*u_objectcolor.rgb, u_objectcolor.a);\n"
		
		//	"    gl_FragDepth=(-(1000.+1.)*(-v_glposition.w)-2.*1000.*1.)/((1000.-1.)*v_glposition.w);\n"
			"    gl_FragDepth=(-(1000.+0.1)*(-v_glposition.w)-2.*1000.*0.1)/((1000.-0.1)*v_glposition.w);\n"
			"}",
			l3d_attr, sizeof(l3d_attr)/sizeof(ShaderVar), l3d_unif, sizeof(l3d_unif)/sizeof(ShaderVar));
		
		ShaderVar bb3d_attr[]=						//Brick Blaster 3D program
		{
			{&GL2_BB3D::a_vertices, "a_vertex", __LINE__},
			{&GL2_BB3D::a_texcoord, "a_texcoord", __LINE__},
			{&GL2_BB3D::a_light, "a_light", __LINE__},
		};
		ShaderVar bb3d_unif[]=
		{
			{&GL2_BB3D::u_matrix, "u_matrix", __LINE__},
			{&GL2_BB3D::u_texture, "u_texture", __LINE__},
		};
		GL2_BB3D::program=LoadShaders(
			"#version 120\n"
			"uniform mat4 u_matrix;\n"
			"attribute vec3 a_vertex;\n"
			"attribute vec2 a_texcoord;\n"
			"attribute float a_light;\n"
			"varying vec2 f_texcoord;\n"
			"varying float v_light;\n"
			"varying vec4 v_fragpos;\n"
			"void main()\n"
			"{\n"
			"    gl_Position=u_matrix*vec4(a_vertex, 1.);\n"
			"    v_fragpos=gl_Position;\n"
			"    gl_Position.z=0.;\n"
			"    f_texcoord=a_texcoord;\n"
			"    v_light=a_light;\n"
			"}",
			"#version 120\n"
			"varying vec2 f_texcoord;\n"
			"varying vec4 v_fragpos;\n"
			"varying float v_light;\n"
			"uniform sampler2D u_texture;\n"
			"void main()\n"
			"{\n"
			"    vec4 color=texture2D(u_texture, f_texcoord);\n"
			"	 gl_FragColor=vec4(color.rgb*v_light, color.a);\n"//alpha is in the texture
		//	"    gl_FragDepth=(-(1000.+1.)*(-v_fragpos.w)-2.*1000.*1.)/((1000.-1.)*v_fragpos.w);\n"
			"    gl_FragDepth=(-(1000.+0.1)*(-v_fragpos.w)-2.*1000.*0.1)/((1000.-0.1)*v_fragpos.w);\n"
			"}",
			bb3d_attr, sizeof(bb3d_attr)/sizeof(ShaderVar), bb3d_unif, sizeof(bb3d_unif)/sizeof(ShaderVar));
		if(!GL2_BB3D::program)
			error(__LINE__);
		
		//glEnable(GL_PROGRAM_POINT_SIZE);
		ShaderVar _3d_attr[]=
		{
			{&GL2_3D::a_vertices, "a_vertex", __LINE__},
		//	{&GL2_3D::a_normals, "a_normal", __LINE__},
			{&GL2_3D::a_texcoord, "a_texcoord", __LINE__}
		};
		ShaderVar _3d_unif[]=
		{
			{&GL2_3D::u_matrix, "u_matrix", __LINE__},
		//	{&GL2_3D::u_pointSize, "u_pointSize", __LINE__},
			{&GL2_3D::u_texture, "u_texture", __LINE__},
			{&GL2_3D::u_isTextured, "u_isTextured", __LINE__},
			{&GL2_3D::u_color, "u_color", __LINE__}
		};
		GL2_3D::program=LoadShaders(						//3D program
			"#version 120\n"
			"uniform mat4 u_matrix;\n"
		//	"uniform float u_pointSize;\n"//can't get u_pointSize, -1
			"attribute vec3 a_vertex;\n"
			"attribute vec2 a_texcoord;\n"
			"varying vec2 f_texcoord;\n"
			"varying vec4 v_fragpos;\n"
			"void main()\n"
			"{\n"
			"    gl_Position=u_matrix*vec4(a_vertex, 1.);\n"
		//	"    gl_Position=u_matrix*vec4(a_vertex, u_pointSize);\n"//got u_pointSize
			"    v_fragpos=gl_Position;\n"
			"    gl_Position.z=0.;\n"
			"    f_texcoord=a_texcoord;\n"
		//	"    gl_PointSize=u_pointSize;\n"
			"    gl_PointSize=20.;\n"
			"}",
			"#version 120\n"
			"varying vec2 f_texcoord;\n"
			"varying vec4 v_fragpos;\n"
			"uniform bool u_isTextured;\n"
			"uniform vec4 u_color;\n"
			"uniform sampler2D u_texture;\n"
			"void main()\n"
			"{\n"
			"    if(u_isTextured)\n"
			"	     gl_FragColor=texture2D(u_texture, f_texcoord);\n"//alpha is in the texture
			"    else\n"
			"        gl_FragColor=u_color;\n"//alpha is in the color
		//	"    gl_FragColor=mix(u_color, texture2D(u_texture, f_texcoord), float(u_isTextured));\n"//samples texture anyway
		
		//	"    gl_FragDepth=(-(1000.+1.)*(-v_fragpos.w)-2.*1000.*1.)/((1000.-1.)*v_fragpos.w);\n"
			"    gl_FragDepth=(-(1000.+0.1)*(-v_fragpos.w)-2.*1000.*0.1)/((1000.-0.1)*v_fragpos.w);\n"
			"}",
			_3d_attr, sizeof(_3d_attr)/sizeof(ShaderVar), _3d_unif, sizeof(_3d_unif)/sizeof(ShaderVar));
		if(!GL2_3D::program)
			error(__LINE__);

		ShaderVar text_attr[]=
		{
			{&GL2_Text::a_coord2d, "coords", __LINE__},
		//	{&GL2_Text::a_coord2d, "coord2d", __LINE__},
		//	{&GL2_Text::a_texcoord, "texcoord", __LINE__}
		};
		ShaderVar text_unif[]=
		{
			{&GL2_Text::u_mytexture, "mytexture", __LINE__},
			{&GL2_Text::u_isTexture, "isTexture", __LINE__},
			{&GL2_Text::u_txtColor, "txtColor", __LINE__},
			{&GL2_Text::u_bkColor, "bkColor", __LINE__}
		};
		GL2_Text::program=LoadShaders(					//Text program
			"#version 120\n"
			"attribute vec4 coords;"
			//"attribute vec2 coord2d;\n"			//coord2d, texcoord
			//"attribute vec2 texcoord;\n"
			"varying vec2 f_texcoord;\n"
			"void main()\n"
			"{\n"
			"    gl_Position=vec4(coords.xy, 0., 1.);\n"
			"    f_texcoord=coords.zw;\n"
			//"    gl_Position=vec4(coord2d, 0., 1.);\n"
			//"    f_texcoord=texcoord;\n"
			"}",
			"#version 120\n"
			"varying vec2 f_texcoord;\n"
			"uniform sampler2D mytexture;\n"	//mytexture, isTexture, txtColor, bkColor
			"uniform vec2 isTexture;\n"
		//	"uniform bool isTexture;\n"
			"uniform vec4 txtColor;\n"
			"uniform vec4 bkColor;\n"
			"void main()\n"
			"{\n"
			"    vec4 region=texture2D(mytexture, f_texcoord);\n"

			"    gl_FragColor=mix((txtColor*region.r+bkColor*region.g), vec4(region.rgb, bkColor.a), isTexture.x);\n"//???fps
		//	"    gl_FragColor=(txtColor*region.r+bkColor*region.g)*isTexture.y+vec4(region.rgb, bkColor.a)*isTexture.x;\n"//38.5fps
		//	"    gl_FragColor=(txtColor*region.r+bkColor*region.g)*float(!isTexture)+vec4(region.rgb, bkColor.a)*float(isTexture);\n"//38.55fps
		//	"    gl_FragColor=txtColor*region.r+bkColor*region.g+region*float(isTexture);\n"//doesn't display texture 41fps, 7.6 fps
		
			//"    if(isTexture.x!=0)\n"
			//"        gl_FragColor=vec4(region.rgb, bkColor.a);\n"//27fps
			//"    else\n"
			//"        gl_FragColor=txtColor*region.r+bkColor*region.g;\n"

			//"    if(isTexture)\n"//texture
			//"        gl_FragColor=vec4(region.rgb, bkColor.a);\n"		//22fps, 3.9~4.1fps
			//"    else if(region.r!=0.)\n"//inside
			//"        gl_FragColor=txtColor;\n"
			//"    else if(region.g!=0.)\n"//background
			//"        gl_FragColor=bkColor;\n"
			//"    else\n"//nodraw
			//"        gl_FragColor=vec4(0.);\n"

			"}",
			text_attr, sizeof(text_attr)/sizeof(ShaderVar), text_unif, sizeof(text_unif)/sizeof(ShaderVar));
		if(!GL2_Text::program)
			error(__LINE__);
		prof_add("Compile shaders");//26.7ms
		
		{
			using namespace resources;
			LoadTexture(sf10, sizeof(sf10), sf10_height, sf10_widths, tx_id_sf10);
			LoadTexture(sf8, sizeof(sf8), sf8_height, sf8_widths, tx_id_sf8);
			LoadTexture(sf7, sizeof(sf7), sf7_height, sf7_widths, tx_id_sf7);
			LoadTexture(sf6, sizeof(sf6), sf6_height, sf6_widths, tx_id_sf6);
		}
		prof_add("Unpack textures");

		{
			float level=1;
			for(int k=15;k>=0;--k)
				light_levels[k]=level, level*=0.8f;//0.8 0.95 0.99
		}
		seed_world(0);
		prof_add("Seed world");
		generate_textures();
		prof_add("Generate textures");
		explore(0, 0, 0, 0);
		prof_add("Generate level");
		gl_setProgram(GL2_Text::program), gl_font_size=gl_setTextSize(4), gl_setTextColor(0xFF000000), gl_setBkColor(0x00FFFFFF);	check(__LINE__);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_CULL_FACE);
		check(__LINE__);
	//	xoshiro128_seed(rand()<<15|rand());
		prof_add("Final preps");
#endif
		return *(int*)&broken;
	}
	void		resize()
	{
		if(h<=0)
			h=1;
		glViewport(0, 0, w, h);
	//	changeFov();
		SetBkMode(ghDC, OPAQUE);
	}
	void		changeFov(){}
	void		clear_screen(){glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);}
	void		show(){SwapBuffers(ghDC);}
	void		finish()
	{
		for(int k=0, kEnd=chunks.size();k<kEnd;++k)//
			glDeleteBuffers(1, &chunks[k].gl_buffer);//
		glDeleteProgram(programID);
		glDeleteBuffers(1, &vbo_cube_vertices);
		glDeleteBuffers(1, &vbo_cube_colors);
		glDeleteBuffers(1, &ibo_cube_elements);
		glDeleteBuffers(1, &vbo_cube_texcoords);
		
		glDeleteProgram(GL2_3D::program);
		glDeleteBuffers(1, &GL2_3D::vertex_buffer);
		
		glDeleteProgram(GL2_Text::program);
		glDeleteBuffers(1, &GL2_Text::buffer);

		wglMakeCurrent(0, 0);
		wglDeleteContext(hRC);
	}
} lgl;
void			switch_engine(int engine)
{
	if(engine!=Engine::mode)
	{
		e->finish();
		if(engine==E_LinearSW)
			e=(Engine*)&lsw;
		else if(engine==E_LinearGL)
			e=(Engine*)&lgl;
		Engine::mode=(EngineType)engine;
		e->initiate();
	}
}
void			switch_engine(){switch_engine(!Engine::mode);}

mat4			*g_mVP=nullptr;
int				collideNo=0;
float			crouch_offset=0;
vec3			collide(vec3 const &p1, vec3 const &d1, vec3 const &p2, vec3 const &d2, char Un, char Dn, char Nn, char Sn, char En, char Wn)
{
	vec3 aa1=p2-d1, aa2=p2+d2;
	const float
		&x1=aa1.x, &y1=aa1.y, &z1=aa1.z,
		&x2=aa2.x, &y2=aa2.y, &z2=aa2.z,
		&x=p1.x, &y=p1.y, &z=p1.z;
	vec3 correction;
	const float maxcorrmagsq=0.05f;//0.125
	if(p1.x>aa1.x&&p1.x<aa2.x&&p1.y>aa1.y&&p1.y<aa2.y&&p1.z>aa1.z&&p1.z<aa2.z)
	{
		float U=z2-z, D=z-z1, N=y2-y, S=y-y1, E=x2-x, W=x-x1;
		float minimum=10;
		enum Face{NONE, UP, DOWN, NORTH, SOUTH, EAST, WEST} face=NONE;
		if(minimum>U&&!Un)
			face=UP, minimum=U;
		if(minimum>D&&!Dn)
			face=DOWN, minimum=D;
		if(minimum>N&&!Nn)
			face=NORTH, minimum=N;
		if(minimum>S&&!Sn)
			face=SOUTH, minimum=S;
		if(minimum>E&&!En)
			face=EAST, minimum=E;
		if(minimum>W&&!Wn)
			face=WEST, minimum=W;
		float C=-0.2f;
		if(debuginfo)
		{
			int linecolor=0xFFFF0000;//
			switch(face)
			{
			case UP:
				label(p2.x+d2.x*0.5f, p2.y+d2.y*0.5f, p2.z+d2.z+0.1f+crouch_offset, "%d: %g", collideNo, z2-z);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, 0.1f, d2.z+0.1f), p2+vec3(0.9f, 0.9f, d2.z+0.1f), linecolor);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, 0.9f, d2.z+0.1f), p2+vec3(0.9f, 0.1f, d2.z+0.1f), linecolor);
				break;
			case DOWN:
				label(p2.x+d2.x*0.5f, p2.y+d2.y*0.5f, p2.z-0.1f+crouch_offset, "%d: %g", collideNo, z1-z);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, 0.1f, -0.1f), p2+vec3(0.9f, 0.9f, -0.1f), linecolor);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, 0.9f, -0.1f), p2+vec3(0.9f, 0.1f, -0.1f), linecolor);
				break;
			case NORTH:
				label(p2.x+d2.x*0.5f, p2.y+d2.y+0.1f, p2.z+d2.z*0.5f+crouch_offset, "%d: %g", collideNo, y2-y);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, d2.y+0.1f, 0.1f), p2+vec3(0.9f, d2.y+0.1f, 0.9f), linecolor);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, d2.y+0.1f, 0.9f), p2+vec3(0.9f, d2.y+0.1f, 0.1f), linecolor);
				break;
			case SOUTH:
				label(p2.x+d2.x*0.5f, p2.y-0.1f, p2.z+d2.z*0.5f+crouch_offset, "%d: %g", collideNo, y1-y);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, -0.1f, 0.1f), p2+vec3(0.9f, -0.1f, 0.9f), linecolor);
				GL2_3D::draw_line(*g_mVP, p2+vec3(0.1f, -0.1f, 0.9f), p2+vec3(0.9f, -0.1f, 0.1f), linecolor);
				break;
			case EAST:
				label(p2.x+d2.x+0.1f, p2.y+d2.y*0.5f, p2.z+d2.z*0.5f+crouch_offset, "%d: %g", collideNo, x2-x);
				GL2_3D::draw_line(*g_mVP, p2+vec3(d2.z+0.1f, 0.1f, 0.1f), p2+vec3(d2.z+0.1f, 0.9f, 0.9f), linecolor);
				GL2_3D::draw_line(*g_mVP, p2+vec3(d2.z+0.1f, 0.1f, 0.9f), p2+vec3(d2.z+0.1f, 0.9f, 0.1f), linecolor);
				break;
			case WEST:
				label(p2.x-0.1f, p2.y+d2.y*0.5f, p2.z+d2.z*0.5f+crouch_offset, "%d: %g", collideNo, x1-x);
				GL2_3D::draw_line(*g_mVP, p2+vec3(-0.1f, 0.1f, 0.1f), p2+vec3(-0.1f, 0.9f, 0.9f), linecolor);
				GL2_3D::draw_line(*g_mVP, p2+vec3(-0.1f, 0.1f, 0.9f), p2+vec3(-0.1f, 0.9f, 0.1f), linecolor);
				break;
			}
		}
		switch(face)
		{
		case UP:
			correction.set(0, 0, z2-z);
			if(correction.mag_sq()>maxcorrmagsq)
				correction.setzero();
			else
			{
				player_v.x*=0.9f, player_v.y*=0.9f;
				if(!kb[VK_SPACE])
					player_v.z*=-0.1f;
				player_v.z=abs(player_v.z);
			}
			return correction;
		case DOWN:	player_v.z*=C; correction.set(0, 0, z1-z);break;
		case NORTH:	player_v.y*=C; correction.set(0, y2-y, 0);break;
		case SOUTH:	player_v.y*=C; correction.set(0, y1-y, 0);break;
		case EAST:	player_v.x*=C; correction.set(x2-x, 0, 0);break;
		case WEST:	player_v.x*=C; correction.set(x1-x, 0, 0);break;
		}
	}
	if(correction.mag_sq()>maxcorrmagsq)
		correction.setzero();
	return correction;
}
float			count_fps()//returns elapsed time in seconds
{
	LARGE_INTEGER li;
	QueryPerformanceFrequency(&li);
	long long freq=li.QuadPart;
	QueryPerformanceCounter(&li);
	static long long nticks=0;
	float delta=float(li.QuadPart-nticks)/freq;
	GUIPrint(ghDC, 0, R.bottom-16, "OpenGL %d.%d, fps=%.10f, T=%.10fms", glMajorVer, glMinorVer, 1/delta, 1000.*delta);
	//linelen=sprintf_s(line, "OpenGL, fps=%.10f, T=%.10fms, dcam=%.10f, fov=%.10f", freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq, dcam, atan(tanfov)*720/_2pi);
	nticks=li.QuadPart;
	return delta;
}
void			count_active_keys()
{
	kp=kb['W']+kb['A']+kb['S']+kb['D']+kb['T']+kb['G']+kb[VK_UP]+kb[VK_DOWN]+kb[VK_LEFT]+kb[VK_RIGHT]
		+kb[VK_ADD]+kb[VK_SUBTRACT]//numpad
		+kb[VK_OEM_PLUS]+kb[VK_OEM_MINUS]
		+kb[VK_RETURN]+kb[VK_BACK];
}
float			delta=0;
int				*rgb2=nullptr, w2=0, h2=0;
vec3			g_corr;
const int
	cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cx+log_cy,
	log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,
	bx=1<<log_bx, by=1<<log_by,//total number of blocks in x & y directions
	dx_mask=level_dx-1, dy_mask=level_dy-1;
void			send_chunk_AO(int ch_idx)
{
	auto &ch=chunks[ch_idx];
	GL2_BB3D::begin_chunk(&ch);
//	ch.offset=GL2_BB3D::bbVertices.size();
	unsigned short gx=ch_idx&dx_mask, gy=ch_idx>>log_level_dx&dy_mask;//chunk idx
	for(auto it=ch.surfaces.begin(), itEnd=ch.surfaces.end();it!=itEnd;++it)
	{
		auto &idx=it->idx;
		int lx=idx&cx_mask, ly=idx>>log_cx&cy_mask,//in-chunk coordinates
			x=ch.x+lx, y=ch.y+ly, z=idx>>log_slice;//world coordinates
		//if(gx==0&&gy==1&&lx==1&&z==125)
		//	int LOL_1=0;
		//if(gx==1&&gy==2&&lx==7&&ly==0&&z==121)
		//	int LOL_1=0;
		bool E=it->DUSNWE&1, W=it->DUSNWE>>1&1, N=it->DUSNWE>>2&1, S=it->DUSNWE>>3&1, U=it->DUSNWE>>4&1, D=it->DUSNWE>>5&1;
		unsigned tx_id=textures[ch.c[idx]%ALL_BLOCKS];
		int
			gyS=ly>0?gy<<log_level_dx:gy>0?(gy-1)<<log_level_dx:-1, gyN=ly<cy_mask?gy<<log_level_dx:gy+1<level_dy?(gy+1)<<log_level_dx:-1, gxW=lx>0?gx:gx>0?gx-1:-1, gxE=lx<cx_mask?gx:gx+1<level_dx?gx+1:-1,//chunk indices, check for -1
			lyS=((ly-1)&cy_mask)<<log_cx, lyN=((ly+1)&cy_mask)<<log_cx, lxW=(lx-1)&cx_mask, lxE=(lx+1)&cx_mask,//always in-bounds
			zD=z>0?(z-1)<<log_slice:-1, zU=z+1<cz?(z+1)<<log_slice:-1,//check for -1
			gxm=gx, gym=gy<<log_level_dx,
			lxm=lx, lym=ly<<log_cx, zm=z<<log_slice;
		byte
			*L_UNW=nullptr, *L_UNm=nullptr, *L_UNE=nullptr,
			*L_UmW=nullptr, *L_Umm=nullptr, *L_UmE=nullptr,
			*L_USW=nullptr, *L_USm=nullptr, *L_USE=nullptr,

			*L_mNW=nullptr, *L_mNm=nullptr, *L_mNE=nullptr,
			*L_mmW=nullptr,					*L_mmE=nullptr,
			*L_mSW=nullptr, *L_mSm=nullptr, *L_mSE=nullptr,

			*L_DNW=nullptr, *L_DNm=nullptr, *L_DNE=nullptr,
			*L_DmW=nullptr, *L_Dmm=nullptr, *L_DmE=nullptr,
			*L_DSW=nullptr, *L_DSm=nullptr, *L_DSE=nullptr;
		Chunk *ch2=nullptr;
		byte skylight=15;
		unsigned short idx2=0;
		if(zU!=-1)
		{
			if(gyN!=-1)
			{
				if(gxW!=-1)	ch2=&chunks[gyN|gxW], idx2=zU|lyN|lxW, L_UNW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
							ch2=&chunks[gyN|gxm], idx2=zU|lyN|lxm, L_UNm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
				if(gxE!=-1)	ch2=&chunks[gyN|gxE], idx2=zU|lyN|lxE, L_UNE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			}
				if(gxW!=-1)	ch2=&chunks[gym|gxW], idx2=zU|lym|lxW, L_UmW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
							ch2=&chunks[gym|gxm], idx2=zU|lym|lxm, L_Umm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
				if(gxE!=-1)	ch2=&chunks[gym|gxE], idx2=zU|lym|lxE, L_UmE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			if(gyS!=-1)
			{
				if(gxW!=-1)	ch2=&chunks[gyS|gxW], idx2=zU|lyS|lxW, L_USW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
							ch2=&chunks[gyS|gxm], idx2=zU|lyS|lxm, L_USm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
				if(gxE!=-1)	ch2=&chunks[gyS|gxE], idx2=zU|lyS|lxE, L_USE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			}
		}
		else
			L_UNW=L_UNm=L_UNE=	L_UmW=L_Umm=L_UmE=	L_USW=L_USm=L_USE=&skylight;

		if(gyN!=-1)
		{
			if(gxW!=-1)	ch2=&chunks[gyN|gxW], idx2=zm|lyN|lxW, L_mNW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
						ch2=&chunks[gyN|gxm], idx2=zm|lyN|lxm, L_mNm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			if(gxE!=-1)	ch2=&chunks[gyN|gxE], idx2=zm|lyN|lxE, L_mNE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
		}
			if(gxW!=-1)	ch2=&chunks[gym|gxW], idx2=zm|lym|lxW, L_mmW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			if(gxE!=-1)	ch2=&chunks[gym|gxE], idx2=zm|lym|lxE, L_mmE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
		if(gyS!=-1)
		{
			if(gxW!=-1)	ch2=&chunks[gyS|gxW], idx2=zm|lyS|lxW, L_mSW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
						ch2=&chunks[gyS|gxm], idx2=zm|lyS|lxm, L_mSm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			if(gxE!=-1)	ch2=&chunks[gyS|gxE], idx2=zm|lyS|lxE, L_mSE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
		}

		if(zD!=-1)
		{
			if(gyN!=-1)
			{
				if(gxW!=-1)	ch2=&chunks[gyN|gxW], idx2=zD|lyN|lxW, L_DNW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
							ch2=&chunks[gyN|gxm], idx2=zD|lyN|lxm, L_DNm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
				if(gxE!=-1)	ch2=&chunks[gyN|gxE], idx2=zD|lyN|lxE, L_DNE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			}
				if(gxW!=-1)	ch2=&chunks[gym|gxW], idx2=zD|lym|lxW, L_DmW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
							ch2=&chunks[gym|gxm], idx2=zD|lym|lxm, L_Dmm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
				if(gxE!=-1)	ch2=&chunks[gym|gxE], idx2=zD|lym|lxE, L_DmE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			if(gyS!=-1)
			{
				if(gxW!=-1)	ch2=&chunks[gyS|gxW], idx2=zD|lyS|lxW, L_DSW=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
							ch2=&chunks[gyS|gxm], idx2=zD|lyS|lxm, L_DSm=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
				if(gxE!=-1)	ch2=&chunks[gyS|gxE], idx2=zD|lyS|lxE, L_DSE=ch2->c[idx2]<=B_GLASS?&ch2->l[idx2]:nullptr;
			}
		}
		if(D)
		{
			float
				V_DSW=vertex_light(L_Dmm, L_DmW, L_DSm, L_DSW),
				V_DSE=vertex_light(L_Dmm, L_DmE, L_DSm, L_DSE),
				V_DNE=vertex_light(L_Dmm, L_DNm, L_DmE, L_DNE),
				V_DNW=vertex_light(L_Dmm, L_DNm, L_DmW, L_DNW);
			GL2_BB3D::push_face_down(x, y, z, 1, 1, tx_id, vec4(V_DSW, V_DSE, V_DNE, V_DNW), V_DSW+V_DNE>V_DSE+V_DNW);
		}
		if(U)
		{
			float
				V_USW=vertex_light(L_Umm, L_UmW, L_USm, L_USW),
				V_USE=vertex_light(L_Umm, L_UmE, L_USm, L_USE),
				V_UNE=vertex_light(L_Umm, L_UNm, L_UmE, L_UNE),
				V_UNW=vertex_light(L_Umm, L_UNm, L_UmW, L_UNW);
			GL2_BB3D::push_face_up(x, y, z+1, 1, 1, tx_id, vec4(V_USW, V_USE, V_UNE, V_UNW), V_UNE+V_USW>V_UNW+V_USE);//correct
		}
		if(S)
		{
			//L_USW L_USm L_USE
			//L_mSW L_mSm L_mSE
			//L_DSW L_DSm L_DSE
			float
				V_DSW=vertex_light(L_mSm, L_mSW, L_DSm, L_DSW),
				V_DSE=vertex_light(L_mSm, L_mSE, L_DSm, L_DSE),
				V_USE=vertex_light(L_mSm, L_mSE, L_USm, L_USE),
				V_USW=vertex_light(L_mSm, L_mSW, L_USm, L_USW);
			GL2_BB3D::push_face_south(x, y, z, 1, 1, tx_id, vec4(V_DSW, V_DSE, V_USE, V_USW), V_DSW+V_USE>V_DSE+V_USW);
		}
		if(N)
		{
			float
				V_DNW=vertex_light(L_mNm, L_mNW, L_DNm, L_DNW),
				V_DNE=vertex_light(L_mNm, L_mNE, L_DNm, L_DNE),
				V_UNE=vertex_light(L_mNm, L_mNE, L_UNm, L_UNE),
				V_UNW=vertex_light(L_mNm, L_mNW, L_UNm, L_UNW);
			GL2_BB3D::push_face_north(x, y+1, z, 1, 1, tx_id, vec4(V_DNW, V_DNE, V_UNE, V_UNW), V_DNW+V_UNE>V_DNE+V_UNW);
		}
		if(W)
		{
			//L_UNW L_UmW L_USW
			//L_mNW L_mmW L_mSW
			//L_DNW L_DmW L_DSW
			float
				V_DSW=vertex_light(L_mmW, L_DmW, L_mSW, L_DSW),
				V_DNW=vertex_light(L_mmW, L_DmW, L_mNW, L_DNW),
				V_UNW=vertex_light(L_mmW, L_UmW, L_mNW, L_UNW),
				V_USW=vertex_light(L_mmW, L_UmW, L_mSW, L_USW);
			GL2_BB3D::push_face_west(x, y, z, 1, 1, tx_id, vec4(V_DSW, V_DNW, V_UNW, V_USW), V_DSW+V_UNW>V_DNW+V_USW);
		}
		if(E)
		{
			float
				V_DSE=vertex_light(L_mmE, L_DmE, L_mSE, L_DSE),
				V_DNE=vertex_light(L_mmE, L_DmE, L_mNE, L_DNE),
				V_UNE=vertex_light(L_mmE, L_UmE, L_mNE, L_UNE),
				V_USE=vertex_light(L_mmE, L_UmE, L_mSE, L_USE);
			GL2_BB3D::push_face_east(x+1, y, z, 1, 1, tx_id, vec4(V_DSE, V_DNE, V_UNE, V_USE), V_DSE+V_UNE>V_DNE+V_USE);
		}
	}
//	ch.count=GL2_BB3D::bbVertices.size()-ch.offset;
	GL2_BB3D::send_chunk();
}
void			send_chunk_blocky(int ch_idx)
{
	auto &ch=chunks[ch_idx];
	GL2_BB3D::begin_chunk(&ch);
	unsigned short gx=ch_idx&dx_mask, gy=ch_idx>>log_level_dx&dy_mask;//chunk idx
	for(auto it=ch.surfaces.begin(), itEnd=ch.surfaces.end();it!=itEnd;++it)
	{
		auto &idx=it->idx;
		int lx=idx&cx_mask, ly=idx>>log_cx&cy_mask,//in-chunk coordinates
			x=ch.x+lx, y=ch.y+ly, z=idx>>log_slice;//world coordinates
		bool E=it->DUSNWE&1, W=it->DUSNWE>>1&1, N=it->DUSNWE>>2&1, S=it->DUSNWE>>3&1, U=it->DUSNWE>>4&1, D=it->DUSNWE>>5&1;
		unsigned tx_id=textures[ch.c[idx]%ALL_BLOCKS];

		int
			gyS=ly>0?gy<<log_level_dx:gy>0?(gy-1)<<log_level_dx:-1, gyN=ly<cy_mask?gy<<log_level_dx:gy+1<level_dy?(gy+1)<<log_level_dx:-1, gxW=lx>0?gx:gx>0?gx-1:-1, gxE=lx<cx_mask?gx:gx+1<level_dx?gx+1:-1,//chunk indices, check for -1
			lyS=((ly-1)&cy_mask)<<log_cx, lyN=((ly+1)&cy_mask)<<log_cx, lxW=(lx-1)&cx_mask, lxE=(lx+1)&cx_mask,//always in-bounds
			zD=z>0?(z-1)<<log_slice:-1, zU=z+1<cz?(z+1)<<log_slice:-1,//check for -1
			gxm=gx, gym=gy<<log_level_dx,
			lxm=lx, lym=ly<<log_cx, zm=z<<log_slice;
		byte *lightinfo;
		float light_level=0;
	//	int ch_idx, idx;
		if(D)
		{
			if(zD!=-1)
			{
				lightinfo=ch.l+(zD|lym|lxm);
				light_level=light_levels[maximum(*lightinfo>>4, *lightinfo&15)];
			}
			else
				light_level=0;
			GL2_BB3D::push_quad(x, y, z, 1, 1, 0, tx_id, true, vec4(light_level));
		}
		if(U)
		{
			if(zU!=-1)
			{
				lightinfo=ch.l+(zU|lym|lxm);
				light_level=light_levels[maximum(*lightinfo>>4, *lightinfo&15)];
			}
			else
				light_level=0;
			GL2_BB3D::push_quad(x, y, z+1, 1, 1, 0, tx_id, false, vec4(light_level));
		}
		if(S)
		{
			if(gyS!=-1)
			{
				lightinfo=chunks[gyS|gxm].l+(zm|lyS|lxm);
				light_level=light_levels[maximum(*lightinfo>>4, *lightinfo&15)];
			}
			else
				light_level=0;
			GL2_BB3D::push_quad(x, y, z, 1, 0, 1, tx_id, false, vec4(light_level));
		}
		if(N)
		{
			if(gyN!=-1)
			{
				lightinfo=chunks[gyN|gxm].l+(zm|lyN|lxm);
				light_level=gyN!=-1?light_levels[maximum(*lightinfo>>4, *lightinfo&15)]:0;
			}
			else
				light_level=0;
			GL2_BB3D::push_quad(x, y+1, z, 1, 0, 1, tx_id, true, vec4(light_level));
		}
		if(W)
		{
			if(gxW!=-1)
			{
				lightinfo=chunks[gym|gxW].l+(zm|lym|lxW);
				light_level=light_levels[maximum(*lightinfo>>4, *lightinfo&15)];
			}
			else
				light_level=0;
			GL2_BB3D::push_quad(x, y, z, 0, 1, 1, tx_id, true, vec4(light_level));
		}
		if(E)
		{
			if(gxE!=-1)
			{
				lightinfo=chunks[gym|gxE].l+(zm|lym|lxE);
				light_level=light_levels[maximum(*lightinfo>>4, *lightinfo&15)];
			}
			else
				light_level=0;
			GL2_BB3D::push_quad(x+1, y, z, 0, 1, 1, tx_id, false, vec4(light_level));
		}
	}
	GL2_BB3D::send_chunk();
}
void			send_chunk_fullbright(int ch_idx)
{
	auto &ch=chunks[ch_idx];
	GL2_BB3D::begin_chunk(&ch);
	unsigned short gx=ch_idx&dx_mask, gy=ch_idx>>log_level_dx&dy_mask;//chunk idx
	for(auto it=ch.surfaces.begin(), itEnd=ch.surfaces.end();it!=itEnd;++it)
	{
		auto &idx=it->idx;
		int lx=idx&cx_mask, ly=idx>>log_cx&cy_mask,//in-chunk coordinates
			x=ch.x+lx, y=ch.y+ly, z=idx>>log_slice;//world coordinates
		bool E=it->DUSNWE&1, W=it->DUSNWE>>1&1, N=it->DUSNWE>>2&1, S=it->DUSNWE>>3&1, U=it->DUSNWE>>4&1, D=it->DUSNWE>>5&1;
		unsigned tx_id=textures[ch.c[idx]%ALL_BLOCKS];
		if(D)
			GL2_BB3D::push_quad(x, y, z, 1, 1, 0, tx_id, true, vec4(1));
		if(U)
			GL2_BB3D::push_quad(x, y, z+1, 1, 1, 0, tx_id, false, vec4(1));
		if(S)
			GL2_BB3D::push_quad(x, y, z, 1, 0, 1, tx_id, false, vec4(1));
		if(N)
			GL2_BB3D::push_quad(x, y+1, z, 1, 0, 1, tx_id, true, vec4(1));
		if(W)
			GL2_BB3D::push_quad(x, y, z, 0, 1, 1, tx_id, true, vec4(1));
		if(E)
			GL2_BB3D::push_quad(x+1, y, z, 0, 1, 1, tx_id, false, vec4(1));
	}
	GL2_BB3D::send_chunk();
}
//inline void extrapolate(
//	float pos_x, float pos_y, float pos_z, float dir_x, float dir_y, float dir_z, float x1, float x2,
//	float *t, float &px1_y, float &px1_z, float &px2_y, float &px2_z)
//{
//	float inv_cx=1/dir_x;
//	t[0]=(x1-pos_x)*inv_cx, px1_y=pos_y+t[0]*dir_y, px1_z=pos_z+t[0]*dir_z;
//	t[1]=(x2-pos_x)*inv_cx, px2_y=pos_y+t[1]*dir_y, px2_z=pos_z+t[1]*dir_z;
//}
//inline void extrapolate(
//	vec3 const &pos, vec3 const &dir, int *p, float x1, float x2,
//	float *t, float &pa_a, float &pa_b, float &pb_a, float &pb_b)
//{
//	float inv_cx=1/dir[p[0]];
//	t[0]=(x1-pos[p[0]])*inv_cx, pa_a=pos[p[1]]+t[0]*dir[p[1]], pa_b=pos[p[2]]+t[0]*dir[p[2]];
//	t[1]=(x2-pos[p[0]])*inv_cx, pb_a=pos[p[1]]+t[1]*dir[p[1]], pb_b=pos[p[2]]+t[1]*dir[p[2]];
//}
inline void		extrapolate(vec3 const &pos, vec3 const &dir, int const *p, float const a1, float const a2, float *t, vec2 &pa1, vec2 &pa2, float const b1, float const b2, float const c1, float const c2)
{
	float inv_cx=1/dir[p[0]];
	t[0]=(a1-pos[p[0]])*inv_cx, pa1.x=pos[p[1]]+t[0]*dir[p[1]], pa1.y=pos[p[2]]+t[0]*dir[p[2]];
	t[1]=(a2-pos[p[0]])*inv_cx, pa2.x=pos[p[1]]+t[1]*dir[p[1]], pa2.y=pos[p[2]]+t[1]*dir[p[2]];
	if(pa1.x<b1||pa1.x>b2||pa1.y<c1||pa1.y>c2)
		t[0]=infinity;
	if(pa2.x<b1||pa2.x>b2||pa2.y<c1||pa2.y>c2)
		t[1]=infinity;
}
//inline int first_positive(float *t, int count)
//{
//	float min_t=t[0];
//	int idx=0;
//	for(int k=1;k<count;++k)
//		if((min_t<0||min_t>t[k])&&t[k]>=0)
//			min_t=t[k], idx=k;
//	return idx;
//}
inline int		first_positive(float *t)
{
	float min_t=-1;
	int idx=-1;
	for(int k=0;k<6;++k)
		if((min_t<0||min_t>t[k])&&t[k]>=0)
			min_t=t[k], idx=k;
	return idx;
}
bool			intersect_ray_AABB(vec3 const &pos, vec3 const &dir, vec3 const &bstart, vec3 const &bend, vec3 &first)
{
	auto &x1=bstart.x, &y1=bstart.y, &z1=bstart.z, &x2=bend.x, &y2=bend.y, &z2=bend.z;
	if(pos.x>=x1&&pos.x<=x2&&pos.y>=y1&&pos.y<=y2&&pos.z>=z1&&pos.z<=z2)
	{
		first=pos;
		return true;
	}
	enum Plane{P_X1, P_X2, P_Y1, P_Y2, P_Z1, P_Z2};
	float t[]={infinity, infinity, infinity, infinity, infinity, infinity};//0:x1	1:x2	2:y1	3:y2	4:z1	5:z2
	vec2 px1_yz, px2_yz, py1_xz, py2_xz, pz1_xy, pz2_xy;
	//float
	//	px1_y, px1_z, px2_y, px2_z,
	//	py1_x, py1_z, py2_x, py2_z,
	//	pz1_x, pz1_y, pz2_x, pz2_y;
	int p_xyz[]={0, 1, 2}, p_yxz[]={1, 0, 2}, p_zxy[]={2, 0, 1};
	if(dir.x)
		extrapolate(pos, dir, p_xyz, x1, x2, t, px1_yz, px2_yz, y1, y2, z1, z2);
	//	extrapolate(pos.x, pos.y, pos.z, dir.x, dir.y, dir.z, x1, x2, t, px1_y, px1_z, px2_y, px2_z);//p1, p2
	else//yz plane
	{
		if(pos.x<x1||pos.x>x2||!dir.y&&!dir.z)
			return false;
		if(dir.y)
			extrapolate(pos, dir, p_yxz, y1, y2, t+2, py1_xz, py2_xz, x1, x2, z1, z2);
		//	extrapolate(pos.y, pos.x, pos.z, dir.y, dir.x, dir.z);
		else//dir.x==0 & dir.y==0
		{
			if(dir.z>0)
			{
				if(pos.z>z2)
					return false;
				if(pos.z>z1)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(pos.x, pos.y, z1);
				return true;
			}
			else
			{
				if(pos.z<z1)
					return false;
				if(pos.z<z2)
				{
					first=pos;//unreachable
					return false;
				}
				first.set(pos.x, pos.y, z2);
				return true;
			}
		}
		if(dir.z)
			extrapolate(pos, dir, p_zxy, z1, z2, t+4, pz1_xy, pz2_xy, x1, x2, y1, y2);
		else//dir.x==0 & dir.z==0
		{
			if(dir.y>0)
			{
				if(pos.y>y2)
					return false;
				if(pos.y>y1)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(pos.x, y1, pos.z);
				return true;
			}
			else
			{
				if(pos.y<y1)
					return nullptr;
				if(pos.y<y2)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(pos.x, y2, pos.z);
				return true;
			}
		}
	//	t[0]=t[1]=infinity;
	}
	//if(dir.y)
	//	extrapolate(pos, dir, p_yxz, y1, y2, t+2, py1_xz, py2_xz);
	//else//horizontal
	//{
	//	if(pos.y<y1||pos.y>y2||!dir.x)
	//		return nullptr;
	//	if(dir.x>0)//looking right
	//	{
	//		if(pos.x>x2)
	//			return nullptr;
	//		if(pos.x>x1)
	//			return pos;//unreachable
	//		return (x1, pos.y);
	//	}
	//	else//looking left
	//	{
	//		if(pos.x<x1)
	//			return nullptr;
	//		if(pos.x<x2)
	//			return pos;
	//		return (x2, pos.y);
	//	}
	//}
	//xyz -> yxz		swap x & y from above
	if(dir.y)
		extrapolate(pos, dir, p_yxz, y1, y2, t+2, py1_xz, py2_xz, x1, x2, z1, z2);
	else//xz plane
	{
		if(pos.y<y1||pos.y>y2||!dir.x&&!dir.z)
			return false;
		if(dir.x)
			extrapolate(pos, dir, p_xyz, x1, x2, t, px1_yz, px2_yz, y1, y2, z1, z2);
		else//dir.x==0 & dir.y==0
		{
			if(dir.z>0)
			{
				if(pos.z>z2)
					return false;
				if(pos.z>z1)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(pos.x, pos.y, z1);
				return true;
			}
			else
			{
				if(pos.z<z1)
					return false;
				if(pos.z<z2)
				{
					first=pos;//unreachable
					return false;
				}
				first.set(pos.x, pos.y, z2);
				return true;
			}
		}
		if(dir.z)
			extrapolate(pos, dir, p_zxy, z1, z2, t+4, pz1_xy, pz2_xy, x1, x2, y1, y2);
		else//dir.y==0 & dir.z==0
		{
			if(dir.x>0)
			{
				if(pos.x>x2)
					return false;
				if(pos.x>x1)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(x1, pos.y, pos.z);
				return true;
			}
			else
			{
				if(pos.x<x1)
					return nullptr;
				if(pos.x<x2)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(x2, pos.y, pos.z);
				return true;
			}
		}
	//	t[2]=t[3]=infinity;
	}
	//xyz -> yxz -> zxy		swap y & z from above
	if(dir.z)
		extrapolate(pos, dir, p_zxy, z1, z2, t+4, pz1_xy, pz2_xy, x1, x2, y1, y2);
	else//xy plane
	{
		if(pos.z<z1||pos.z>z2||!dir.x&&!dir.y)
			return false;
		if(dir.x)
			extrapolate(pos, dir, p_xyz, x1, x2, t, px1_yz, px2_yz, y1, y2, z1, z2);
		else//dir.x==0 & dir.z==0
		{
			if(dir.y>0)
			{
				if(pos.y>y2)
					return false;
				if(pos.y>y1)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(pos.x, y1, pos.z);
				return true;
			}
			else
			{
				if(pos.y<y1)
					return false;
				if(pos.y<y2)
				{
					first=pos;//unreachable
					return false;
				}
				first.set(pos.x, y2, pos.z);
				return true;
			}
		}
		if(dir.y)
			extrapolate(pos, dir, p_yxz, y1, y2, t+2, py1_xz, py2_xz, x1, x2, z1, z2);
		else//dir.y==0 & dir.z==0
		{
			if(dir.x>0)
			{
				if(pos.x>x2)
					return false;
				if(pos.x>x1)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(x1, pos.y, pos.z);
				return true;
			}
			else
			{
				if(pos.x<x1)
					return nullptr;
				if(pos.x<x2)
				{
					first=pos;//unreachable
					return true;
				}
				first.set(x2, pos.y, pos.z);
				return true;
			}
		}
	//	t[4]=t[5]=infinity;
	}
	switch(first_positive(t))
	{
	case P_X1:first.set(x1, px1_yz.x, px1_yz.y);return true;
	case P_X2:first.set(x2, px2_yz.x, px2_yz.y);return true;
	case P_Y1:first.set(py1_xz.x, y1, py1_xz.y);return true;
	case P_Y2:first.set(py2_xz.x, y2, py2_xz.y);return true;
	case P_Z1:first.set(pz1_xy.x, pz1_xy.y, z1);return true;
	case P_Z2:first.set(pz2_xy.x, pz2_xy.y, z2);return true;
	}
	return false;
}
byte*			containing_block(vec3 const &pos, vec3 const &dir, int &wx, int &wy, int &z)
{
	float X0=float(bx>>1), Y0=float(by>>1);
	vec3 start;
	if(intersect_ray_AABB(pos, dir, vec3(-X0, -Y0, 0), vec3(X0-1, Y0-1, (float)cz-1), start))
	{
		wx=(int)floor(start.x), wy=(int)floor(start.y), z=(int)floor(start.z);//
		return get(wx, wy, z);
	}
	return nullptr;
}
struct			Voxel
{
	int wx, wy, z;
	byte *p;
	Voxel(int wx, int wy, int z, byte *p):wx(wx), wy(wy), z(z), p(p){}
};
#if 0
bool update_alg=true;//
vec3 g_pos, g_dir;//
std::vector<Voxel> voxels;//
typedef std::pair<float, int> TC;//t, color
std::vector<TC> g_t;//
byte*			cast_ray_debug(vec3 const &pos, vec3 const &dir, int &wx, int &wy, int &z)
{
	voxels.clear(), g_t.clear();
	float X0=float(bx>>1), Y0=float(by>>1);
	vec3 start;
	if(intersect_ray_AABB(pos, dir, vec3(-X0, -Y0, 0), vec3(X0-1, Y0-1, (float)cz-1), start))
	{
		int dx=(dir.x>0)-(dir.x<0), dy=(dir.y>0)-(dir.y<0), dz=(dir.z>0)-(dir.z<0);//what if the ray is aligned to grid (one of them is zero)
		wx=(int)floor(start.x), wy=(int)floor(start.y), z=(int)floor(start.z);
		if(wx<-X0||wx>=X0||wy<-Y0||wy>=Y0||z<0||z>=cz)
			return nullptr;
		byte *current=get(wx, wy, z);
		voxels.push_back(Voxel(wx, wy, z, current));//
		if(current&&*current>B_AIR)
			return current;
	//	vec3 xy, yz, xz;
	//	int next_x=wx+dx, next_y=wy+dy, next_z=z+dz;//https://github.com/francisengelmann/fast_voxel_traversal/blob/master/main.cpp
		vec3										//http://www.cse.yorku.ca/~amana/research/grid.pdf
			t(
				dir.x?((dx>0?ceil(start.x):floor(start.x))-pos.x)/dir.x:infinity,
				dir.y?((dy>0?ceil(start.y):floor(start.y))-pos.y)/dir.y:infinity,
				dir.z?((dz>0?ceil(start.z):floor(start.z))-pos.z)/dir.z:infinity),
		//	t(dir.x?(next_x-pos.x)/dir.x:infinity, dir.y?(next_y-pos.y)/dir.y:infinity, dir.z?(next_z-pos.z)/dir.z:infinity),
			dt(dir.x?dx/dir.x:infinity, dir.y?dy/dir.y:infinity, dir.z?dz/dir.z:infinity),
			diff;
		bool neg_ray=false;
		while(wx>=-X0&&wx<X0&&wy>=-Y0&&wy<Y0&&z>=0&&z<cz)
		{
			if(t.x<t.y)
			{
				if(t.x<t.z)
				{
					wx+=dx, t.x+=dt.x;
					g_t.push_back(TC(t.x, 0xFF000000|rand()<<15|rand()));
				}
				else
				{
					z+=dz, t.z+=dt.z;
					g_t.push_back(TC(t.z, 0xFF000000|rand()<<15|rand()));
				}
			}
			else
			{
				if(t.y<t.z)
				{
					wy+=dy, t.y+=dt.y;
					g_t.push_back(TC(t.y, 0xFF000000|rand()<<15|rand()));
				}
				else
				{
					z+=dz, t.z+=dt.z;
					g_t.push_back(TC(t.z, 0xFF000000|rand()<<15|rand()));
				}
			}
			current=get(wx, wy, z);
			voxels.push_back(Voxel(wx, wy, z, current));//
			if(current&&*current>B_AIR)
				return current;
		}
	}
	return nullptr;
}
#endif
byte*			cast_ray(vec3 const &pos, vec3 const &dir, int &wx, int &wy, int &z, int &wx0, int &wy0, int &z0)//http://www.cse.yorku.ca/~amana/research/grid.pdf
{
	float X0=float(bx>>1), Y0=float(by>>1);
	vec3 start;
	if(intersect_ray_AABB(pos, dir, vec3(-X0, -Y0, 0), vec3(X0-1, Y0-1, (float)cz-1), start))
	{
		int dx=(dir.x>0)-(dir.x<0), dy=(dir.y>0)-(dir.y<0), dz=(dir.z>0)-(dir.z<0);//what if the ray is aligned to grid (one of them is zero)
		wx0=wx=(int)floor(start.x), wy0=wy=(int)floor(start.y), z0=z=(int)floor(start.z);
		if(wx<-X0||wx>=X0||wy<-Y0||wy>=Y0||z<0||z>=cz)
			return nullptr;
		byte *current=get(wx, wy, z);
		if(current&&*current>B_AIR)
			return current;
		vec3 t( dir.x?((dx>0?ceil(start.x):floor(start.x))-pos.x)/dir.x:infinity,
				dir.y?((dy>0?ceil(start.y):floor(start.y))-pos.y)/dir.y:infinity,
				dir.z?((dz>0?ceil(start.z):floor(start.z))-pos.z)/dir.z:infinity),
			dt(dir.x?dx/dir.x:infinity, dir.y?dy/dir.y:infinity, dir.z?dz/dir.z:infinity);
		while(wx>=-X0&&wx<X0&&wy>=-Y0&&wy<Y0&&z>=0&&z<cz)
		{
			wx0=wx, wy0=wy, z0=z;
			if(t.x<t.y)
			{
				if(t.x<t.z)
					wx+=dx, t.x+=dt.x;
				else
					z+=dz, t.z+=dt.z;
			}
			else
			{
				if(t.y<t.z)
					wy+=dy, t.y+=dt.y;
				else
					z+=dz, t.z+=dt.z;
			}
			current=get(wx, wy, z);
			if(current&&*current>B_AIR)
				return current;
		}
	}
	return nullptr;
}
void			render()
{
	prof_add("Render entry");
	LARGE_INTEGER li;
	QueryPerformanceFrequency(&li);
	long long freq=li.QuadPart;
	QueryPerformanceCounter(&li);
	long long t1=li.QuadPart, t2=__rdtsc();
	float time_ms=(float)(t1-t_start)/freq*1000;
	const float max_delta=0.02f;
	if(delta>max_delta)
		delta=max_delta;
	prof_add("Calculate time");
	
	bool crouching=kb[VK_CONTROL];
	if(noclip)
	{
		if(timer)
		{
				 if(kb[VK_SHIFT]){	 if(kb['W'])	cam.moveFastForward();
									 if(kb['A'])	cam.moveFastLeft();
									 if(kb['S'])	cam.moveFastBack();
									 if(kb['D'])	cam.moveFastRight();
									 if(kb['T'])	cam.moveFastUp();
									 if(kb['G'])	cam.moveFastDown();}
			else{					 if(kb['W'])	cam.moveForward();
									 if(kb['A'])	cam.moveLeft();
									 if(kb['S'])	cam.moveBack();
									 if(kb['D'])	cam.moveRight();
									 if(kb['T'])	cam.moveUp();
									 if(kb['G'])	cam.moveDown();}
				 if(kb[VK_UP])						cam.turnUp();
				 if(kb[VK_DOWN])					cam.turnDown();
				 if(kb[VK_LEFT])					cam.turnLeft();
				 if(kb[VK_RIGHT])					cam.turnRight();
				 if(kb[VK_ADD		]||kb[VK_RETURN	]||kb[VK_OEM_PLUS	])	cam.zoomIn();
				 if(kb[VK_SUBTRACT	]||kb[VK_BACK	]||kb[VK_OEM_MINUS	])	cam.zoomOut();
		}
	}
	else//player movement
	{
		if(timer)
		{
			if(kb[VK_SHIFT]&&!crouching)
			{
				if(kb['W'])		cam.playerRunForward(player_v, 50*delta);
				if(kb['A'])		cam.playerRunLeft(player_v, 50*delta);
				if(kb['S'])		cam.playerRunBack(player_v, 50*delta);
				if(kb['D'])		cam.playerRunRight(player_v, 50*delta);
			}
			else
			{
				if(kb['W'])		cam.playerWalkForward(player_v, 50*delta);
				if(kb['A'])		cam.playerWalkLeft(player_v, 50*delta);
				if(kb['S'])		cam.playerWalkBack(player_v, 50*delta);
				if(kb['D'])		cam.playerWalkRight(player_v, 50*delta);
			}
			if(kb[VK_UP])		cam.turnUp();
			if(kb[VK_DOWN])		cam.turnDown();
			if(kb[VK_LEFT])		cam.turnLeft();
			if(kb[VK_RIGHT])	cam.turnRight();
			if(kb[VK_ADD		]||kb[VK_RETURN	]||kb[VK_OEM_PLUS	])	cam.zoomIn();
			if(kb[VK_SUBTRACT	]||kb[VK_BACK	]||kb[VK_OEM_MINUS	])	cam.zoomOut();
		}
	}
	crouch_offset=crouching*player_d.z*0.5f;
	prof_add("Controls");

	e->clear_screen();
	mat4 mView=cam.viewMatrix(crouch_offset),
	//	mProj=perspective(cam.tanfov, float(w)/h, 1, 1000),
		mProj=perspective(cam.tanfov, float(w)/h, 0.1f, 1000),
		mVP=mProj*mView;
	g_mVP=&mVP;
	prof_add("Matrix");

	vec3 p0=cam.p;
	p0.z-=crouch_offset;
	if(!noclip)
	{
		player_v.z-=9.8f*delta;//add acceleration due to gravity
		cam.p+=player_v*delta;//update player position
	}
	prof_add("Simulate");
	float
		x1, x2, dx,
		y1, y2, dy,
		z1, z2, dz;
	std::vector<vec3> corrections;//
	if(!noclip||noclipCollisionsOn)
	{
		vec3 p1=cam.p-player_d*0.5f, p2;
		p1.z+=player_d.z*0.5f-player_eye_height;
		p2=p1+player_d;
		p2.z-=crouch_offset;
		//const float
		//	&x1=p1.x, &x2=p2.x, &dx=player_d.x,
		//	&y1=p1.y, &y2=p2.y, &dy=player_d.y,
		//	&z1=p1.z, &z2=p2.z, &dz=player_d.z;
		x1=p1.x, x2=p2.x, dx=player_d.x,
		y1=p1.y, y2=p2.y, dy=player_d.y,
		z1=p1.z, z2=p2.z, dz=player_d.z-crouch_offset;
		int fx1=(int)floor(x1), fx2=(int)floor(x2), fy1=(int)floor(y1), fy2=(int)floor(y2), fz1=(int)floor(z1), fz2=(int)floor(z2);
		int nbh_x=fx2+1-fx1, nbh_y=fy2+1-fy1, nbh_z=fz2+1-fz1;
		bool old_correction=true;//
		collideNo=0;
		for(int kz=fz1;kz<=fz2;++kz)
			for(int ky=fy1;ky<=fy2;++ky)
				for(int kx=fx1;kx<=fx2;++kx)
				{
					auto bk=get(kx, ky, kz);
					if(bk&&*bk>B_AIR)
					{
						byte *Un=get(kx, ky, kz+1), *Dn=get(kx, ky, kz-1), *Nn=get(kx, ky+1, kz), *Sn=get(kx, ky-1, kz), *En=get(kx+1, ky, kz), *Wn=get(kx-1, ky, kz);
						vec3 correction=collide(p1, vec3(dx, dy, dz), vec3((float)kx, (float)ky, (float)kz), vec3(1), Un&&*Un>B_AIR, Dn&&*Dn>B_AIR, Nn&&*Nn>B_AIR, Sn&&*Sn>B_AIR, En&&*En>B_AIR, Wn&&*Wn>B_AIR);
						corrections.push_back(correction);//
						cam.p+=correction, p1+=correction;
						if(correction.mag_sq())
						{
							if(old_correction)
								g_corr=correction, old_correction=false;
							else
								g_corr+=correction;
						}
						++collideNo;
					}
				}
	}
	mView=cam.viewMatrix(crouch_offset), mVP=mProj*mView;
	prof_add("Collision detection");
	if(debuginfo)
	{
		vec3 p1=cam.p-player_d*0.5f, p2;
		p1.z+=player_d.z*0.5f-player_eye_height;
		p2=p1+player_d;
		p2.z-=crouch_offset;
		x1=p1.x, x2=p2.x, dx=player_d.x,
		y1=p1.y, y2=p2.y, dy=player_d.y,
		z1=p1.z, z2=p2.z, dz=player_d.z-crouch_offset;
		GL2_3D::draw_AABB(mVP, x1, x2, y1, y2, z1, z2, 0xFFFF00FF);//draw player bounding box
		//int bound_color=0xFFFF00FF;
		//GL2_3D::draw_line(mVP, vec3(x1, y1, z1), vec3(x2, y1, z1), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x2, y1, z1), vec3(x2, y2, z1), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x2, y2, z1), vec3(x1, y2, z1), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x1, y2, z1), vec3(x1, y1, z1), bound_color);//
		//
		//GL2_3D::draw_line(mVP, vec3(x1, y1, z1), vec3(x1, y1, z2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x2, y1, z1), vec3(x2, y1, z2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x2, y2, z1), vec3(x2, y2, z2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x1, y2, z1), vec3(x1, y2, z2), bound_color);//
		//
		//GL2_3D::draw_line(mVP, vec3(x1, y1, z2), vec3(x2, y1, z2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x2, y1, z2), vec3(x2, y2, z2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x2, y2, z2), vec3(x1, y2, z2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(x1, y2, z2), vec3(x1, y1, z2), bound_color);//
		
		//bound_color=0xFF000000;
		//float mx1, mx2, my1, my2, mz1, mz2;
		//GL2_3D::draw_line(mVP, vec3(mx1, my1, mz1), vec3(mx2, my1, mz1), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx2, my1, mz1), vec3(mx2, my2, mz1), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx2, my2, mz1), vec3(mx1, my2, mz1), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx1, my2, mz1), vec3(mx1, my1, mz1), bound_color);//
		//
		//GL2_3D::draw_line(mVP, vec3(mx1, my1, mz1), vec3(mx1, my1, mz2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx2, my1, mz1), vec3(mx2, my1, mz2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx2, my2, mz1), vec3(mx2, my2, mz2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx1, my2, mz1), vec3(mx1, my2, mz2), bound_color);//
		//
		//GL2_3D::draw_line(mVP, vec3(mx1, my1, mz2), vec3(mx2, my1, mz2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx2, my1, mz2), vec3(mx2, my2, mz2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx2, my2, mz2), vec3(mx1, my2, mz2), bound_color);//
		//GL2_3D::draw_line(mVP, vec3(mx1, my2, mz2), vec3(mx1, my1, mz2), bound_color);//
	}
	GL2_3D::draw_AABB(mVP, -(float)draw_distance, (float)draw_distance, -(float)draw_distance, (float)draw_distance, 0, (float)cz, 0xFF000000);//draw map bounding box

	//3D Example
#if 0
	static float angle=0;
	angle+=0.2f*!paused;
	mat4
		anim=
			rotate(mat4(1), angle*3.0f*torad, vec3(1, 0, 0)) *  // X axis
			rotate(mat4(1), angle*2.0f*torad, vec3(0, 1, 0)) *  // Y axis
			rotate(mat4(1), angle*4.0f*torad, vec3(0, 0, 1)),   // Z axis
	//	anim=rotate(mat4(1), angle*torad, vec3(0, 1, 0)),
		model=translate(mat4(1), vec3(0, 0, -2)),//-4
		view=matrixFPSViewRH(cam, a.y, a.x-_pi),
	//	view=lookAt(cam, cam+vec3(cax*cay, -sax*cay, say), vec3(0, 1, 0)),//X
	//	view=lookAt(vec3(0, 2, 0), vec3(0, 0, -4), vec3(0, 1, 0)),
		projection=perspective(tanfov, (float)w/h, 0.1f, 1000),
	//	projection=perspective(fov*torad*w/h, (float)w/h, 0.1f, 1000),//BROKEN
	//	projection=perspective(fov*torad, (float)w/h, 0.1f, 1000),
		mvp=projection*view*model*anim;
	
	gl_setProgram(programID);										check(__LINE__);
	glUniformMatrix4fv(uniform_mvp, 1, GL_FALSE, &mvp[0][0]);		check(__LINE__);

	select_texture(tx_id_sf10, uniform_texture);
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, tx_id_sf10);
	//glUniform1i(uniform_texture, 0);//GL_TEXTURE

	glUniform1f(uniform_fade, 1);
	//glUniform1f(uniform_fade, 0.5f*(1+sin(time_ms*0.01f)));
	glEnableVertexAttribArray(attribute_coord3d);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_vertices);
	glVertexAttribPointer(
		attribute_coord3d,	// attribute
		3,					// number of elements per vertex, here (x,y,z)
		GL_FLOAT,			// the type of each element
		GL_FALSE,			// take our values as-is
		0,					// no extra data between each position
		0					// offset of first element
	);
	glEnableVertexAttribArray(attribute_texcoord);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_texcoords);
	glVertexAttribPointer(
		attribute_texcoord, // attribute
		2,                  // number of elements per vertex, here (x,y)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_cube_elements);
	int idx_buf_size;
	glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &idx_buf_size);
	glDrawElements(GL_TRIANGLES, idx_buf_size>>1, GL_UNSIGNED_SHORT, 0);
	//glDrawArrays(GL_TRIANGLES, 0, 3);
	glDisableVertexAttribArray(attribute_coord3d);
	glDisableVertexAttribArray(attribute_texcoord);
#endif
	//light test
#if 0
	gl_setProgram(GL2_L3D::program);
	vec3 lightpos(10, 10, 10), spherepos(-10, -10, -10);
	mat4
		mView=matrixFPSViewRH(cam, a.y, a.x-_pi),
		mProj=perspective(tanfov, float(w)/h, 0.1f, 1000.f),
		vp=mProj*mView;
	{//draw the sphere
		mat4 model=translate(mat4(1), spherepos);
		mat3 m_normal=normalMatrix(model);
		mat4 mvp=vp*model;
		glUniformMatrix4fv(GL2_L3D::u_vpmatrix, 1, GL_FALSE, mvp.data());			check(__LINE__);
		glUniformMatrix4fv(GL2_L3D::u_modelmatrix, 1, GL_FALSE, model.data());		check(__LINE__);
		glUniformMatrix3fv(GL2_L3D::u_normalmatrix, 1, GL_FALSE, m_normal.data());	check(__LINE__);
		send_color(GL2_L3D::u_objectcolor, 0xFF0000FF);								check(__LINE__);
		send_color_rgb(GL2_L3D::u_lightcolor, 0xFFFFFF);							check(__LINE__);
		glUniform3fv(GL2_L3D::u_lightpos, 1, &lightpos[0]);							check(__LINE__);
		glUniform3fv(GL2_L3D::u_viewpos, 1, &cam[0]);								check(__LINE__);
	//	glBindVertexArray(s_VAO);													check(__LINE__);
		glPointSize(10);															check(__LINE__);
		
		glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);												check(__LINE__);
		glVertexAttribPointer(GL2_L3D::a_vertices, 3, GL_FLOAT, GL_FALSE, 6<<2, 0);				check(__LINE__);
		glEnableVertexAttribArray(GL2_L3D::a_vertices);											check(__LINE__);
		glVertexAttribPointer(GL2_L3D::a_normals, 3, GL_FLOAT, GL_FALSE, 6<<2, (void*)(3<<2));	check(__LINE__);
		glEnableVertexAttribArray(GL2_L3D::a_normals);											check(__LINE__);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereEBO);										check(__LINE__);

		glDrawElements(GL_TRIANGLES, s_vcount, GL_UNSIGNED_INT, 0);	check(__LINE__);
	//	glDrawArrays(GL_POINTS, 0, s_vcount);						check(__LINE__);
	}
#endif
	//3D test
#if 0
	static bool not_sent=true;
	if(not_sent)
	{
		not_sent=false;

		w2=128, h2=128;
		int size=w2*h2;
		rgb2=new int[size];//allocated once
		for(int k=0;k<size;++k)
			rgb2[k]=0x7F000000|(rand()<<15|rand())&0x00FFFFFF;//alpha in texture

		GL2_3D::begin();

		//int black=0xFF000000;
		vec3 zero(0, 0, 0);
		GL2_3D::push_line_segment(zero, vec3(10, 0, 0), 0xFF0000FF);//x red
		GL2_3D::push_line_segment(zero, vec3(0, 10, 0), 0xFF00FF00);//y green
		GL2_3D::push_line_segment(zero, vec3(0, 0, 10), 0xFFFF0000);//z blue
		GL2_3D::push_triangle(vec3(-0.5f, 0.25f-0.5f, 3), vec3(0.5f, 0.25f-0.5f, 3), vec3(0,  0.25f+0.5f, 3), 0xFF000FF0);

		GL2_3D::begin_transparent();
		for(int z=-10;z<0;z+=2)
			GL2_3D::push_square(-10, 10, -10, 10, (float)z, rgb2, w2, h2);
	//	GL2_3D::push_square(-10, 10, -10, 10, -20, rgb2, w2, h2);
	//	GL2_3D::push_square(-10, 10, -10, 10, 20, rgb2, w2, h2);
		GL2_3D::push_triangle(vec3(-0.5f, -0.5f, 2), vec3(0.5f, -0.5f, 2), vec3(0,  0.5f, 2), 0x7F0000FF);
		GL2_3D::push_triangle(vec3(-0.5f, -0.5f, 4), vec3(0.5f, -0.5f, 4), vec3(0,  0.5f, 4), 0x7F00FF00);
		GL2_3D::push_triangle(vec3(-0.5f, -0.5f, 8), vec3(0.5f, -0.5f, 8), vec3(0,  0.5f, 8), 0x7FFF0000);
		GL2_3D::end();//sends buffers to GPU
	}
	GL2_3D::draw();
#endif
	//2D test, region test
#if 0
//	GL2_2D::set_region(0, w, 0, h);
	GL2_2D::set_region(w>>2, 3*(w>>2), h>>2, 3*(h>>2));
	GL2_2D::use_region();
//	glViewport(w>>2, h>>2, w>>1, h>>1);
	{
		double *curve=new double[w];
		double
			VX=0, DX=20, Xstart=VX-DX*0.5, Xsample=DX/w,
			VY=0, DY=DX*h/w, Yend=VY+DY*0.5, inv_Ysample=h/DY;
		for(int k=0;k<w;++k)
			curve[k]=cos(Xstart+Xsample*k);

		GL2_2D::pen_color=0xFFFF0000;
		GL2_2D::curve_start();
		for(int k=0;k<w;++k)
			GL2_2D::curve_point((float)k, float((Yend-curve[k])*inv_Ysample));
		//	GL2_2D::curve_point((float)k, (float)k);

		gl_setProgram(GL2_2D::program);				check(__LINE__);
		send_color(GL2_2D::u_color, 0xFFFF00FF);	check(__LINE__);
		//for(int k=0;k<100;++k)
		//	GL2_2D::draw_line(rand()%w, rand()%h, rand()%w, rand()%h);
		GL2_2D::draw_line(10, 10, w-10, h-10);
		GL2_2D::draw_line(w-10, 10, 10, h-10);
		send_color(GL2_2D::u_color, 0x7F7F7F7F);
		GL2_2D::draw_rectangle(w>>1, w*3>>2, h>>1, h*3>>2);
		GL2_2D::draw_rectangle_hollow(w>>2, w>>1, h>>2, h>>1);
		delete[] curve;
	}
	GL2_2D::draw_curve();
	GL2_2D::drop_region();
//	glViewport(0, 0, w, h);
#endif
	//brick blaster
#if 1
	static bool update=false;
	static int ch_idx=0;
	if(!lights_done)
	{
		int idx1=start_gy<<log_level_dx|start_gx;
		calculate_lights_v3(false);
		if(lights_done)
			update=true, ch_idx=0;
		//if(lights_done)
		//	int LOL_1=0;
		int idx2=start_gy<<log_level_dx|start_gx;
		for(int k=idx1;k<idx2;++k)
			send_chunk_AO(k);
	//	not_sent=true;
		prof_add("Calculate lights v3");
	}
	if(update)
	{
		if(ch_idx<(int)chunks.size())
			send_chunk_AO(ch_idx++);
		else
			update=false;
	}
	//if(lighting&&light_iterations<15)
	//{
	//	calculate_lights_iterate();
	//	++light_iterations;
	//	not_sent=true;
	//	prof_add("Iterate lights");
	//}
	if(not_sent)
	{
		not_sent=false;
	//	GL2_BB3D::begin();
		GL2_3D::begin();
		if(lighting==2)//ambient occlusion
		{
			for(int k=0, kEnd=chunks.size();k<kEnd;++k)
			{
				std::vector<int> LOL_1(1, 0);
				send_chunk_AO(k);
			}
		}
		else if(lighting==1)//blocky lighting
		{
			for(int k=0, kEnd=chunks.size();k<kEnd;++k)
				send_chunk_blocky(k);
		}
		else//fullbright
		{
			for(int k=0, kEnd=chunks.size();k<kEnd;++k)
				send_chunk_fullbright(k);
		}
		vec3 zero(0, 0, 0);
		GL2_3D::push_line_segment(zero, vec3(1000, 0, 0), 0xFF0000FF);//x red
		GL2_3D::push_line_segment(zero, vec3(0, 1000, 0), 0xFF00FF00);//y green
		GL2_3D::push_line_segment(zero, vec3(0, 0, 1000), 0xFFFF0000);//z blue
		GL2_3D::begin_transparent();
	//	GL2_BB3D::end();
		GL2_3D::end();
	}
	GL2_BB3D::draw(mVP);//TODO: don't draw chunks outside view
	GL2_3D::draw(mVP);
#endif
	{
		int wx, wy, z, wx0, wy0, z0;
		vec3 p=cam.p, dir=cam.viewray();
		p.z-=crouch_offset;
		//if(update_alg)
		//{
		//	g_pos=cam.p, g_dir=dir;
			byte *selection=cast_ray(p, dir, wx, wy, z, wx0, wy0, z0);
		//}
	/*	if(voxels.size())
		{
			int selection_color=0xFF7F7F7F;
			float a=1.02f, b=0.51f;
			for(int k=0, kEnd=voxels.size()-1;k<kEnd;++k)//
			{
				auto &vk=voxels[k];
				vec3 p((float)vk.wx, (float)vk.wy, (float)vk.z);
			//	vec3 p((float)vk.wx-0.01f, (float)vk.wy-0.01f, (float)vk.z-0.01f);
				//GL2_3D::draw_line(mVP, p, p+vec3(a, a, a), selection_color);
				//GL2_3D::draw_line(mVP, p+vec3(a, 0, 0), p+vec3(0, a, a), selection_color);
				//GL2_3D::draw_line(mVP, p+vec3(a, a, 0), p+vec3(0, 0, a), selection_color);
				//GL2_3D::draw_line(mVP, p+vec3(0, a, 0), p+vec3(a, 0, a), selection_color);
				GL2_3D::draw_AABB(mVP, p.x, p.x+1, p.y, p.y+1, p.z, p.z+1, selection_color);
			//	GL2_3D::draw_AABB(mVP, p.x, p.x+a, p.y, p.y+a, p.z, p.z+a, selection_color);
				label(p.x+b, p.y+b, p.z+b, "(%d,%d, %d)", vk.wx, vk.wy, vk.z);
			}
			selection_color=0xFF000000;
			auto &vk2=*voxels.rbegin();
			vec3 p2((float)vk2.wx-0.01f, (float)vk2.wy-0.01f, (float)vk2.z-0.01f);
			GL2_3D::draw_AABB(mVP, p2.x, p2.x+a, p2.y, p2.y+a, p2.z, p2.z+a, selection_color);
			{
				vec3 lp1, lp2=g_pos;
				for(int k=0, kEnd=g_t.size();k<kEnd;++k)
				{
					lp1=lp2, lp2=g_pos+g_t[k].first*g_dir;
				//	lp1.z++, lp2.z+=k+1;
					GL2_3D::draw_line(mVP, lp1, lp2, g_t[k].second);
					label(lp2.x, lp2.y, lp2.z, "%g", g_t[k].first);
				}
			}
		}//*/
		if(selection)
		{
			int selection_color=0xFF000000;//0xFF7F7F7F
			vec3 p((float)wx-0.01f, (float)wy-0.01f, (float)z-0.01f);
			float a=1.02f;
			//glPointSize(2);
			GL2_3D::draw_AABB(mVP, p.x, p.x+a, p.y, p.y+a, p.z, p.z+a, selection_color);
			//GL2_3D::draw_line(mVP, p,				p+vec3(a, 0, 0),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(a, 0, 0), p+vec3(a, a, 0),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(a, a, 0), p+vec3(0, a, 0),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(0, a, 0), p,					selection_color);
			//
			//GL2_3D::draw_line(mVP, p+vec3(0, 0, 0), p+vec3(0, 0, a),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(a, 0, 0), p+vec3(a, 0, a),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(a, a, 0), p+vec3(a, a, a),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(0, a, 0), p+vec3(0, a, a),	selection_color);
			//
			//GL2_3D::draw_line(mVP, p+vec3(0, 0, a), p+vec3(a, 0, a),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(a, 0, a), p+vec3(a, a, a),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(a, a, a), p+vec3(0, a, a),	selection_color);
			//GL2_3D::draw_line(mVP, p+vec3(0, a, a), p+vec3(0, 0, a),	selection_color);
			//glPointSize(1);
		}//*/
	}
	prof_add("Render world");
	//Text test
#if 1
	glDisable(GL_DEPTH_TEST);
	if(debuginfo)
	{
		vec3 p=p0+vec3(cam.cax*cam.cay, cam.sax*cam.cay, cam.say);//draw axes at crosshair
		float a=0.1f;
		GL2_3D::draw_line(mVP, p, p+vec3(a, 0, 0), 0xFF0000FF);//x red
		GL2_3D::draw_line(mVP, p, p+vec3(0, a, 0), 0xFF00FF00);//y green
		GL2_3D::draw_line(mVP, p, p+vec3(0, 0, a), 0xFFFF0000);//z blue
		label(p.x+a, p.y, p.z+crouch_offset, "East");
		label(p.x, p.y+a, p.z+crouch_offset, "North");
		label(p.x, p.y, p.z+a+crouch_offset, "Up");
	}
	prof_add("Axes & bounding box");
	gl_setProgram(GL2_Text::program);		check(__LINE__);
	//display texture test
#if 0
	int
		w2=20, h2=10,
	//	w2=w, h2=h,
	//	w2=w>>1, h2=h>>1,
		size2=w2*h2, *rgb2=(int*)malloc(w2*h2<<2);
	for(int k=0;k<size2;++k)
	//	rgb2[k]=unsigned char(k/w2*255/h2)<<16|unsigned char((k%w2)*255/w2);
		rgb2[k]=0x00FFFFFF&-((k/w2&1)^(k%w2&1));
	//	rgb2[k]=k;
	//{
	//	auto p=(unsigned char*)(rgb2+k);
	//	p[0]=p[1]=p[2]=rand();
	//}
	//	rgb2[k]=rand()<<15|rand();
	//rgb2[0]=0xFFFF00FF;
	display_texture(0, w>>1, 0, h>>1, rgb2, w2, h2, 0x50);
	//display_texture(10, w-10, 20, h-20, rgb2, w2, h2, 0x50);
	free(rgb2);
#endif
	//print close block coordinates
#if 1
	{
		const int
			cx_mask=cx-1, cy_mask=cy-1, log_slice=log_cx+log_cy,
			log_bx=log_level_dx+log_cx, log_by=log_level_dy+log_cy,
			bx=1<<log_bx, by=1<<log_by,//total number of blocks in x & y directions
			dx_mask=level_dx-1, dy_mask=level_dy-1;
		for(int k=0, kEnd=chunks.size();k<kEnd;++k)
		{
			auto &ch=chunks[k];
			unsigned short gx=k&dx_mask, gy=k>>log_level_dx&dy_mask;//chunk idx
			for(auto it=ch.surfaces.begin(), itEnd=ch.surfaces.end();it!=itEnd;++it)
			{
				auto &idx=it->idx;
				int lx=idx&cx_mask, ly=idx>>log_cx&cy_mask,//in-chunk coordinates
					x=ch.x+lx, y=ch.y+ly, z=idx>>log_slice;//world coordinates
				vec3 p((float)x, (float)y, (float)z);
				if((p-cam.p).mag_sq()<25)//print idx within 5 blocks
					label(x+0.5f, y+0.5f, z+0.5f, "(%d, %d), (%d, %d, %d)", gx, gy, lx, ly, z);//
				//	label(x+0.5f, y+0.5f, z+0.5f, "%01X %04X", k, idx);//
			}
		}
	}
	prof_add("Block coordinates");
#endif
	text_show();
	int fontH=gl_getFontHeight();
	if(debuginfo)
		for(int k=0, kEnd=corrections.size();k<kEnd;++k)//
		{
			auto &ck=corrections[k];
			gl_print(0, k*fontH, 0, "(%f, %f, %f)", ck.x, ck.y, ck.z);//not visible
		}
	gl_setTextColor(txtColor);		check(__LINE__);
	gl_setBkColor(bkColor);			check(__LINE__);
	//for(int k=0;k<1000;++k)
	//	gl_print(rand()%w, rand()%h, "Stress test. Sample text.");
	//gl_print(10, 20, 0, "Hello World,\tpress \'Q\'");
	if(debuginfo)
		gl_print(w>>1, h>>1, 0, "lighting=%d", lighting);
	//for(int k=1;k<20;++k)//tab origin test
	//	gl_print(10, 20+gl_getFontHeight()*k, 10*k, "Hello World,\tpress \'Q\'");
	//gl_print(rand()%w, rand()%h, 0, "ABCDEFGH");
	//gl_print(10, 20, 0, "Hello World");
	//gl_print(0, 20+gl_getFontHeight(), 0, "2222222\t2");

	if(debuginfo)//DEBUG
	{
		gl_print(0, h-16-11*fontH, 0, "light: %s, chunk(%d, %d), column(%d, %d), update=%d", lights_done?"DONE":"NOT DONE", start_gx, start_gy, start_lx, start_ly, (int)update);
	//	gl_print(0, h-16-10*fontH, 0, "light_iterations=%d", (int)light_iterations);
		gl_print(0, h-16-8*fontH, 0, "p.z=%g", cam.p.z-player_eye_height);
		gl_print(0, h-16-6*fontH, 0, "v.x=%g", player_v.x);
		gl_print(0, h-16-5*fontH, 0, "v.y=%g", player_v.y);
		gl_print(0, h-16-4*fontH, 0, "v.z=%g", player_v.z);
	}
	else//draw crosshair
	{
		int X0=w>>1,Y0=h>>1, clen=20, cwid=2;//10
		static bool oneshot=true;
		if(oneshot)
			GL2_2D::drop_region(), oneshot=false;
		//GL2_2D::set_region(0, w, 0, h);
		//GL2_2D::use_region();
		GL2_2D::set_color(0x7FFFFFFF);//0x7F7F7F7F
		//gl_setProgram(GL2_2D::program);
		//send_color(GL2_2D::u_color, 0x7F7F7F7F);
		//GL2_2D::draw_rectangle(0, w, 0, h);//
		GL2_2D::draw_rectangle(X0-clen, X0+clen, Y0-cwid, Y0+cwid);
		GL2_2D::draw_rectangle(X0-cwid, X0+cwid, Y0-clen, Y0-cwid);
		GL2_2D::draw_rectangle(X0-cwid, X0+cwid, Y0+cwid, Y0+clen);
		cwid=1;
		GL2_2D::set_color(0x7F000000);
		GL2_2D::draw_rectangle(X0-clen, X0+clen, Y0-cwid, Y0+cwid);
		GL2_2D::draw_rectangle(X0-cwid, X0+cwid, Y0-clen, Y0-cwid);
		GL2_2D::draw_rectangle(X0-cwid, X0+cwid, Y0+cwid, Y0+clen);
		//GL2_2D::draw_line(X0-10, Y0-1, X0+10, Y0-1);
		//GL2_2D::draw_line(X0-10, Y0  , X0+10, Y0  );
		//GL2_2D::draw_line(X0-10, Y0+1, X0+10, Y0+1);
		//GL2_2D::draw_line(X0-1, Y0-10, X0-1, Y0+10);
		//GL2_2D::draw_line(X0  , Y0-10, X0  , Y0+10);
		//GL2_2D::draw_line(X0+1, Y0-10, X0+1, Y0+10);
	}
	if(help)
	{
		gl_setTextColor(0xFF000000);		check(__LINE__);
		gl_setBkColor(0xFFFFFFFF);			check(__LINE__);
		const char *msg[]=
		{
			"W/A/S/D/T/G: move",
			"Arrows: turn",
			"Space: jump",
			"Mouse1: remove block",
			"Esc: toggle mouse control",
		//	"Mouse1: toggle mouse control",
			"Wheel/+/-/Enter/Backspace: change FOV",
			"Ctrl wheel: change text size",
			"Shift wheel: change noclip speed",
			"N: noclip",
			"F: change light type",
			"H: recalculate lights",
			"R: reset position & FOV",
			"Q: toggle timer",//
			"F1: print help",
			"F2: toggle noclip collisions",
			"F3: toggle debug info",
			"F5: new world",
		};
		int hx=w>>2, hy=h>>2;
		for(int k=0, kEnd=sizeof(msg)>>2;k<kEnd;++k)
			gl_print(hx, hy+k*fontH, hx, msg[k]);
		gl_setTextColor(txtColor);			check(__LINE__);
		gl_setBkColor(bkColor);				check(__LINE__);
	}
	gl_print(0, h-16-3*fontH, 0, "seed = %016llX", world_seed);
//	gl_print(0, h-16-3*fontH, 0, "v(%g, %g, %g)", player_v.x, player_v.y, player_v.z);
	gl_print(0, h-16-2*fontH, 0, "fov=%g, dcam=%g, pos(%g, %g, %g/%g), angles(%g, %g)", 2*todeg*atan(cam.tanfov), cam.dcam, cam.p.x, cam.p.y, cam.p.z-crouch_offset, cam.p.z-player_eye_height, cam.a.x, cam.a.y);
	QueryPerformanceCounter(&li);
	static long long nticks=0;
	gl_print(0, h-16-fontH, 0, "OpenGL %s, fps=%.10f, T=%.10fms", GLversion, freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq);
	nticks=li.QuadPart;
	check(__LINE__);
	glEnable(GL_DEPTH_TEST);
#endif
	//long long t_start=__rdtsc();
	//Sleep(100);
	//int dt=int(__rdtsc()-t_start);

	e->show();
	
	//GUIPrint(ghDC, w>>1, h>>1, dt);
	
	if(broken)
	{
		int code=*(int*)&broken, line=((int*)&broken)[1];
		if(line)
		{
			if(code)
				GUIPrint(ghDC, w>>1, h>>1, "ERROR: code %d, line %d", code, line);
			else
				GUIPrint(ghDC, w>>1, h>>1, "ERROR: line %d", line);
		}
		else
			GUIPrint(ghDC, w>>1, h>>1, "ERROR: code %d", code);
	}
	prof_print();
	delta=count_fps();
}
void			__stdcall TimerCallback(unsigned timerID, unsigned msg, unsigned long user, unsigned long _1, unsigned long _2)
{
	count_fps();
}
void			__stdcall TimerProc(HWND hWindow, unsigned msg, unsigned timerID, unsigned long elapsed)
{
	render();
	//count_fps();
}
long			__stdcall WndProc(HWND__ *hWnd, unsigned message, unsigned wParam, long lParam)
{
	switch(message)
	{
	//case WM_CREATE:
	//	break;
	case WM_PAINT:
		GetClientRect(hWnd, &R);
		if(h!=R.bottom-R.top||w!=R.right-R.left)
		{
			h=R.bottom-R.top, w=R.right-R.left, centerP.x=X0=w/2, centerP.y=Y0=h/2;
			ClientToScreen(hWnd, &centerP);
			e->resize();
		}
		if(!timer)
			render();
		break;
	case WM_EXITSIZEMOVE:
		centerP.x=w>>1, centerP.y=h>>1;
		ClientToScreen(ghWnd, &centerP);
		return 0;
	case WM_TIMER:
		render();
	//	count_active_keys();
		if(!kp)
			KillTimer(hWnd, 0), timer=false;
		break;
	case WM_LBUTTONDOWN:
		if(!drag)//enter mouse control
		{
			drag=true;
			ShowCursor(false);
			mouseP0.x=short(lParam), mouseP0.y=short(lParam>>16);
			ClientToScreen(hWnd, &mouseP0);
			SetCursorPos(centerP.x, centerP.y);
		}
		else//remove block
		{
			vec3 pos=cam.p, dir=cam.viewray();
			pos.z-=crouch_offset;
			int wx, wy, z, ch_idx, gx, gy, b_idx, lx, ly, wx0, wy0, z0;
			byte *selection=cast_ray(pos, dir, wx, wy, z, wx0, wy0, z0);
			if(selection&&(ch_idx=get(wx, wy, z, gx, gy, b_idx, lx, ly))!=-1&&*selection!=B_BEDROCK)
			{
				if(*selection==B_TORCH)
				{
				}
				else
				{
					byte *block=selection, *lightinfo=chunks[ch_idx].l+b_idx,
						*bE=nullptr, *bW=nullptr, *bN=nullptr, *bS=nullptr, *bU=nullptr, *bD=nullptr,
						*lE=nullptr, *lW=nullptr, *lN=nullptr, *lS=nullptr, *lU=nullptr, *lD=nullptr;
					get_world(wx+1, wy, z, bE, lE);
					get_world(wx-1, wy, z, bW, lW);
					get_world(wx, wy+1, z, bN, lN);
					get_world(wx, wy-1, z, bS, lS);
					get_world(wx, wy, z+1, bU, lU);
					get_world(wx, wy, z-1, bD, lD);
					char
						skylight=maximum(
							bE&&*bE<=B_GLASS?(*lE&15)-1:0, bW&&*bW<=B_GLASS?(*lW&15)-1:0,
							bN&&*bN<=B_GLASS?(*lN&15)-1:0, bS&&*bS<=B_GLASS?(*lS&15)-1:0,
							bU&&*bU<=B_GLASS?*lU&15:0, bD&&*bD<=B_GLASS?(*lD&15)-1:0),
						torchlight=maximum(
							bE&&*bE<=B_GLASS?*lE>>4:0, bW&&*bW<=B_GLASS?*lW>>4:0,
							bN&&*bN<=B_GLASS?*lN>>4:0, bS&&*bS<=B_GLASS?*lS>>4:0,
							bU&&*bU<=B_GLASS?*lU>>4:0, bD&&*bD<=B_GLASS?*lD>>4:0)-1;
					*lightinfo=torchlight<<4|skylight;
				}
				*selection=B_AIR;//remove the block
				update_lights(ch_idx, b_idx);

				update_surfaces(gx, gy);
				send_chunk_AO(ch_idx);
			//	chunks[ch_idx].update_surfaces(gx+1>=level_dx?nullptr:&chunks[gy<<log_level_dx|(gx+1)], gx-1<0?nullptr:&chunks[gy<<log_level_dx|(gx-1)], gy+1>=level_dy?nullptr:&chunks[(gy+1)<<log_level_dx|gx], gy-1<0?nullptr:&chunks[(gy-1)<<log_level_dx|gx]);
				if(lx==0)//also update neighbor chunks if on chunk edge
				{
					if(gx-1>0)
					{
						update_surfaces(gx-1, gy);
						send_chunk_AO(gy<<log_level_dx|(gx-1));
					}
				}
				else if(lx==cx_mask)
				{
					if(gx+1>0)
					{
						update_surfaces(gx+1, gy);
						send_chunk_AO(gy<<log_level_dx|(gx+1));
					}
				}
				if(ly==0)
				{
					if(gy-1>0)
					{
						update_surfaces(gx, gy-1);
						send_chunk_AO((gy-1)<<log_level_dx|gx);
					}
				}
				else if(ly==cy_mask)
				{
					if(gy+1>0)
					{
						update_surfaces(gx, gy+1);
						send_chunk_AO((gy+1)<<log_level_dx|gx);
					}
				}
			}
			if(!timer)
				render();
		}
		break;
	//	return 0;
	case WM_RBUTTONDOWN://place block
		{
			vec3 pos=cam.p, dir=cam.viewray();
			pos.z-=crouch_offset;
			int wx, wy, z, ch_idx, gx, gy, b_idx, lx, ly, wx0, wy0, z0;
			byte *selection=cast_ray(pos, dir, wx, wy, z, wx0, wy0, z0);
			if(selection&&(ch_idx=get(wx0, wy0, z0, gx, gy, b_idx, lx, ly))!=-1)
			{
				auto &ch=chunks[ch_idx];
				ch.c[b_idx]=B_DIRT, ch.l[b_idx]=0;//place the block
			//	update_lights(ch_idx, b_idx);

				update_surfaces(gx, gy);
				send_chunk_AO(ch_idx);
			//	chunks[ch_idx].update_surfaces(gx+1>=level_dx?nullptr:&chunks[gy<<log_level_dx|(gx+1)], gx-1<0?nullptr:&chunks[gy<<log_level_dx|(gx-1)], gy+1>=level_dy?nullptr:&chunks[(gy+1)<<log_level_dx|gx], gy-1<0?nullptr:&chunks[(gy-1)<<log_level_dx|gx]);
				if(lx==0)//also update neighbor chunks if on chunk edge
				{
					if(gx-1>0)
					{
						update_surfaces(gx-1, gy);
						send_chunk_AO(gy<<log_level_dx|(gx-1));
					}
				}
				else if(lx==cx_mask)
				{
					if(gx+1>0)
					{
						update_surfaces(gx+1, gy);
						send_chunk_AO(gy<<log_level_dx|(gx+1));
					}
				}
				if(ly==0)
				{
					if(gy-1>0)
					{
						update_surfaces(gx, gy-1);
						send_chunk_AO((gy-1)<<log_level_dx|gx);
					}
				}
				else if(ly==cy_mask)
				{
					if(gy+1>0)
					{
						update_surfaces(gx, gy+1);
						send_chunk_AO((gy+1)<<log_level_dx|gx);
					}
				}
			}
			if(!timer)
				render();
		}
		break;
	case WM_MOUSEMOVE://task manager sometimes causes WM_MOUSEMOVE as it updates
		if(drag&&(bypass=!bypass))
		{
			cam.turnMouse(lParam);
			SetCursorPos(centerP.x, centerP.y);
			if(!timer)
				render();
		}
		if(!drag)
			mouseP0.x=short(lParam), mouseP0.y=short(lParam>>16);
		break;
	case WM_MOUSEWHEEL:
		{
			char mw_forward=short(wParam>>16)>0;
			if(kb[VK_SHIFT])//shift wheel		change cam speed
			{
					 if(mw_forward)	cam.faster();
				else				cam.slower();
			}
			//	dcam*=mw_forward*2+!mw_forward*0.5f;
			else if(kb[VK_CONTROL])//ctrl wheel		change text size
				gl_font_size=gl_setTextSize(gl_font_size+mw_forward-!mw_forward);
			else//wheel
			{
					 if(mw_forward)	cam.zoomIn();
				else				cam.zoomOut();
			}
			if(!timer)
				render();
		}
		break;
	case WM_KEYDOWN:case WM_SYSKEYDOWN:
		switch(wParam)
		{
		case 'W':case 'A':case 'S':case 'D':case 'T':case 'G':case VK_UP:case VK_DOWN:case VK_LEFT:case VK_RIGHT:
		case VK_ADD:case VK_SUBTRACT:case VK_OEM_PLUS:case VK_OEM_MINUS:case VK_RETURN:case VK_BACK:
			if(!kb[wParam])
				kb[wParam]=true, count_active_keys();
			if(!timer)
				SetTimer(hWnd, 0, USER_TIMER_MINIMUM, 0), timer=true;
			break;

	/*	case 'W':
			if(!timer)//DEBUG
			{
				cam.moveForward();
			//	cam.playerWalkForward(player_v, 0.01);
				render();//
			}
			break;
		case 'A':
			if(!timer)//
			{
				cam.moveLeft();
				render();//
			}
			break;
		case 'S':
			if(!timer)//
			{
				cam.moveBack();
				render();//
			}
			break;
		case 'D':
			if(!timer)//
			{
				cam.moveRight();
				render();//
			}
			break;//*/

		case VK_SPACE://jump
			if(!noclip)
			{
				player_v.z+=5;
				if(!timer)
					render();
			}
			break;
		case VK_CONTROL://crouch
			if(!noclip)
			{
				if(!timer)
					render();
			}
			break;
		//case '1':
		//	update_alg=!update_alg;
		//	if(!timer)
		//		render();
		//	break;
		case 'E':
			if(!timer)
				render();
			break;
		case VK_ESCAPE:
			ShowCursor(drag);
			if(drag=!drag)
			{
			//	mouseP0.x=short(lParam), mouseP0.y=short(lParam>>16);
				ClientToScreen(hWnd, &mouseP0);
				SetCursorPos(centerP.x, centerP.y);
			}
			else
				SetCursorPos(mouseP0.x, mouseP0.y);
			break;
		case VK_F1:
			help=!help;
			if(!timer)
				render();
			break;
		case VK_F2:
			noclipCollisionsOn=!noclipCollisionsOn;
			if(!timer)
				render();
			break;
		case VK_F3:
			debuginfo=!debuginfo;
			if(!timer)
				render();
			break;
		case VK_F5://new world
			{
				seed_world(__rdtsc());
				explore(0, 0, 0, 0);
				lights_done=false;
			//	not_sent=true;
				if(!timer)
					render();
			}
			break;
		case 'F'://change lighting type
			//if(!lighting)//was fullbright
			//	calculate_lights_v3(true);
			//	calculate_lights_rain();
			//	calculate_lights();
			lighting=(lighting+1)%3;
			not_sent=true;
			if(!timer)
				render();
			break;
		//case 'L'://iterate light
		//	if(lighting&&light_iterations<15)
		//	{
		//		calculate_lights_iterate();
		//		++light_iterations;
		//		not_sent=true;
		//		if(!timer)
		//			render();
		//	}
		//	break;
		//case 'K'://calculate all lights v1
		//	if(lighting)
		//	{
		//		calculate_lights();
		//		prof_add("Calculate lights v1");//223ms
		//		not_sent=true;
		//		if(!timer)
		//			render();
		//	}
		//	break;
		//case 'J'://calculate all lights v2
		//	if(lighting)
		//	{
		//		calculate_lights_v2();
		//		prof_add("Calculate lights v2");//208ms
		//		not_sent=true;
		//		if(!timer)
		//			render();
		//	}
		//	break;
		case 'H'://calculate lights iteratively v3
			if(lighting)
			{
				lights_done=false;
				calculate_lights_v3(true);
			//	not_sent=true;
				if(!timer)
					render();
			}
			break;
		case 'N'://toggle noclip
			noclip=!noclip;
			if(noclip)
				player_v.setzero();
			//paused=!paused;
			//if(timer=!timer)
			//	SetTimer(hWnd, 0, USER_TIMER_MINIMUM, TimerProc);
			//else
			//	KillTimer(hWnd, 0);
			if(!timer)
				render();
			break;
		case 'R':
			cam.reset();
		//	GL2_3D::reset_cam();
			if(!timer)
				render();
			break;
		//case VK_RETURN:
		//	gl_font_size=gl_setTextSize(gl_font_size+1);
		//	break;
		//case VK_BACK:
		//	gl_font_size=gl_setTextSize(gl_font_size-1);
		//	break;
		case 'Q':
			paused=!paused;
			//if(timer=!timer)//100 fps
			//	timerID=timeSetEvent(10, 10, TimerCallback, 0, TIME_PERIODIC|TIME_KILL_SYNCHRONOUS);
			//else
			//	timeKillEvent(timerID);
			if(timer=!timer)//64 fps
				SetTimer(hWnd, 0, USER_TIMER_MINIMUM, TimerProc);
			else
				KillTimer(hWnd, 0);
			//for(;;)//~350 fps
			//{
			//	render();
			//	count_fps();
			//}
			break;
		}
		kb[wParam]=true;
		break;
	case WM_KEYUP:case WM_SYSKEYUP:
		kb[wParam]=false;
		switch(wParam)
		{
		case 'W':case 'A':case 'S':case 'D':case 'T':case 'G':case VK_UP:case VK_DOWN:case VK_LEFT:case VK_RIGHT:
		case VK_ADD:case VK_SUBTRACT:case VK_OEM_PLUS:case VK_OEM_MINUS:case VK_RETURN:case VK_BACK:
			count_active_keys();
			break;
		}
		break;
	case WM_CLOSE:
		PostQuitMessage(0);
		break;
	}
	return DefWindowProcA(hWnd, message, wParam, lParam);
}
int				__stdcall WinMain(HINSTANCE__ *hInstance, HINSTANCE__*, char*, int nCmdShow)
{
	tagWNDCLASSEXA wndClassEx={sizeof(tagWNDCLASSEXA), CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, WndProc, 0, 0, hInstance, LoadIconA(0, (char*)0x00007F00), LoadCursorA(0, (char*)0x00007F00), (HBRUSH__*)(COLOR_WINDOW+1), 0, "New format", 0};
	RegisterClassExA(&wndClassEx);
	ghWnd=CreateWindowExA(0, wndClassEx.lpszClassName, "Brick Blaster 20200513", WS_CAPTION|WS_SYSMENU|WS_THICKFRAME|WS_MINIMIZEBOX|WS_MAXIMIZEBOX|WS_CLIPCHILDREN, CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, 0, 0, hInstance, 0);
	ShowWindow(ghWnd, nCmdShow);

	//	prof_start();
		path.resize(MAX_PATH);
		GetModuleFileNameA(0, &path[0], MAX_PATH);
		for(int k=path.size()-1;k>=0;--k)
		{
			if(path[k]=='/'||path[k]=='\\')
			{
				path.resize(k+1);
				break;
			}
		}
	//	prof_add("path");

		LARGE_INTEGER li;
		QueryPerformanceCounter(&li);
		t_start=li.QuadPart;

		ghDC=GetDC(ghWnd);
		e=(Engine*)&lgl;
		GetClientRect(ghWnd, &R);
		if(h!=R.bottom-R.top||w!=R.right-R.left)
		{
			h=R.bottom-R.top, w=R.right-R.left, centerP.x=X0=w/2, centerP.y=Y0=h/2;
			ClientToScreen(ghWnd, &centerP);
			e->initiate();
		}
	//	prof_add("init");

	tagMSG msg;
	for(;GetMessageA(&msg, 0, 0, 0);)TranslateMessage(&msg), DispatchMessageA(&msg);

		e->finish();
		ReleaseDC(ghWnd, ghDC);

	return msg.wParam;
}