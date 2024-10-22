// Raytracing_n.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include "pch.h"
#include "common.h"
#include <ctime>
#include <chrono>
#include <thread>
#include <mutex>
#include "sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include "stdlib.h"
#include "material.h"
#include "moving_sphere.h"
#include "aarect.h"
#include "triangle.h"
#include "box.h"
#include "constant_medium.h"
#include "bvh.h"
#include "pdf.h"
#include "model.h"
#include "brdf.h"
#include "env_sphere.h"
#include "teapot.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define RaysBackgroundY

using namespace std;

const int thread_count = 8;
int **colors;
int finishedPixel = 0;
int curfinishedcount;
int threadCount = 1;
int threadfinishedcount = 0;
int nx = 1000;
int ny = 1000;
int ns = 50;
int maxDepth = 50;
int sceneid = 2;		// 0: bunny 1: teapot 2: balls 3: orennayar balls 4: jade bunny 5: finals 6: soldier 7: flatnormal_bunny
mutex g_lock;
ofstream outfile("..\\results\\20210709_balls.ppm", ios_base::out);

inline vec3 de_nan(const vec3& c) {
	vec3 temp = c;
	if (!(temp[0] == temp[0])) temp[0] = 0;
	if (!(temp[1] == temp[1])) temp[1] = 0;
	if (!(temp[2] == temp[2])) temp[2] = 0;
	return temp;
}
int bad_sample = 0;
vec3 color(const ray& r, hitable *world, hitable *light_shape, int *depth)
{
	hit_record hrec;
	if (world->hit(r, 0.001, numeric_limits<float>::max(), hrec))
	{
		scatter_record srec;
		vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
		float pdf_val = 0;
		if (*depth < maxDepth && hrec.mat_ptr->scatter(r, hrec, srec))
		//if (hrec.mat_ptr->scatter(r, hrec, srec))
		{
			if (srec.is_specular)
			{
				*depth += 1;
				return srec.attenuation * color(srec.specular_ray, world, light_shape, depth);
			}
			else
			{
				ray scattered;
				pdf_val = 0;
				if (((hitable_list*)light_shape)->list_size > 0)
				{
					hitable_pdf plight(light_shape, hrec.p);
					mixture_pdf p(&plight, srec.pdf_ptr); 
					while (pdf_val == 0)
					{
						scattered = ray(hrec.p, p.generate(r.direction()), r.time());
						pdf_val = p.value(r.direction(), scattered.direction());
					}
				}
				else
				{
					scattered = ray(hrec.p, srec.pdf_ptr->generate(r.direction()), r.time());
					pdf_val = srec.pdf_ptr->value(r.direction(), scattered.direction());
				}
				

				delete srec.pdf_ptr;
				*depth += 1;
				return emitted + srec.attenuation * hrec.mat_ptr->scattering_pdf(r, hrec, scattered) * color(scattered, world, light_shape, depth) / pdf_val;
				//
			}
		}else
		{
			return emitted;
		}
	}
	else
	{
		return vec3(0.0);
	}
}

void random_scene(hitable **scene, camera **cam, hitable **hlist, float aspect) {

	vec3 lookfrom(-10, 6, -15);
	vec3 lookat(0, 0, 0);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int n = 500;
	hitable **list = new hitable*[n + 1];
	texture *checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9, 0.9, 0.9)));
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checker));
	int i = 1;
	for (int a = -11; a < 11; a++)
	{
		for (int b = -11; b < 11; b++)
		{
			float choose_mat = drand48();
			vec3 center(a + 0.9*drand48(), 0.2, b + drand48());
			if ((center - vec3(4, 0.2, 0)).length() > 0.9)
			{
				if (choose_mat < 0.8) {
					//list[i++] = new sphere(center, 0.2, new lambertian(vec3(drand48()*drand48(), drand48()*drand48(), drand48()*drand48())));
					list[i++] = new moving_sphere(center, center + vec3(0, 0.5 * drand48(), 0), 0.0, 1.0, 0.2, new lambertian(new constant_texture(vec3(0.5, 0.5, 0.5))));
				}
				else if (choose_mat < 0.95)
				{
					list[i++] = new sphere(center, 0.2,
						new metal(vec3(0.5*(1 + drand48()), 0.5*(1 + drand48()), 0.5*(1 + drand48())), 0.5*drand48()));
				}
				else {
					list[i++] = new sphere(center, 0.2, new dielectric(1.5));
				}
			}
		}
	}
	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
	list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

	// environment 
	int tx, ty, tn;
	hitable** env_shape = new hitable*[6];
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_1\\Front.jpg", &tx, &ty, &tn, 0);							// front
	env_shape[0] = new xy_rect(-100, 100, -100, 100, -100, new diffuse_light(new image_texture(tex_data, tx, ty)));
	tex_data = stbi_load("..\\contents\\environment_map\\sky_1\\Back.jpg", &tx, &ty, &tn, 0);						// front
	env_shape[1] = new flip_normals(new xy_rect(-100, 100, -100, 100, 100, new diffuse_light(new image_texture(tex_data, tx, ty))));		// back
	tex_data = stbi_load("..\\contents\\environment_map\\sky_1\\Left.jpg", &tx, &ty, &tn, 0);						// front
	env_shape[2] = new yz_rect(-100, 100, -100, 100, 100, new diffuse_light(new image_texture(tex_data, tx, ty)));						// left
	tex_data = stbi_load("..\\contents\\environment_map\\sky_1\\Right.jpg", &tx, &ty, &tn, 0);						// front
	env_shape[3] = new flip_normals(new yz_rect(-100, 100, -100, 100, -100, new diffuse_light(new image_texture(tex_data, tx, ty))));		// right
	tex_data = stbi_load("..\\contents\\environment_map\\sky_1\\Top.jpg", &tx, &ty, &tn, 0);						// front
	env_shape[4] = new flip_normals(new xz_rect(-100, 100, -100, 100, 100, new diffuse_light(new image_texture(tex_data, tx, ty))));		// top
	tex_data = stbi_load("..\\contents\\environment_map\\sky_1\\Bottom.jpg", &tx, &ty, &tn, 0);						// front
	env_shape[5] = new xz_rect(-100, 100, -100, 100, -100, new diffuse_light(new image_texture(tex_data, tx, ty)));						// bottom

	for (size_t j = 0; j < 6; j++)
	{
		list[i++] = env_shape[j];
	}
	//env_shape[6] = new sphere(vec3(0, 6, 0), 0.1, new diffuse_light(new constant_texture(vec3(100, 100, 100))));

	*hlist = new hitable_list(env_shape, 6);
	//delete tex_data;

	*scene = new hitable_list(list, i);
	//return new hitable_list(list, i);
}

hitable *two_spheres() {
	texture *checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9, 0.9, 0.9)));
	int n = 50;
	hitable ** list = new hitable*[n + 1];
	list[0] = new sphere(vec3(0, -10, 0), 10, new lambertian(checker));
	list[1] = new sphere(vec3(0, 10, 0), 10, new lambertian(checker));

	return new hitable_list(list, 2);
}

hitable *two_perlin_spheres() {
	texture* pertext = new noise_texture(1);
	hitable **list = new hitable*[2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	return new hitable_list(list, 2);
}

hitable *earth_shpere() {
	int nx, ny, nn;
	unsigned char* tex_data = stbi_load(".\\Contents\\earthmap.jpg", &nx, &ny, &nn, 0);
	material* mat = new diffuse_light(new image_texture(tex_data, nx, ny));
	hitable **list = new hitable*[2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.9, 0.9, 0.9))));
	list[1] = new sphere(vec3(0, 2, 0), 2, mat);
	return new hitable_list(list, 2);
}

hitable *simple_light() {
	texture *pertext = new noise_texture(4);
	hitable **list = new hitable*[3];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	//list[2] = new sphere(vec3(0, 7, 0), 2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
	list[2] = new xy_rect(3, 5, 1, 3, -2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
	return new hitable_list(list, 3);
}

void cornell_box(hitable **scene, camera **cam, hitable **hlist, float aspect) {

	vec3 lookfrom(300, 500, -800);
	vec3 lookat(300, 278, 200);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	brdf* aluminum_brdf = new brdf();
	hitable **list = new hitable*[200];
	material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
	material *red_200 = new lambertian(new constant_texture(vec3(0.93, 0.6, 0.6)));
	material *white = new lambertian(new constant_texture(vec3(1)));
	material *green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
	material *blue = new lambertian(new constant_texture(vec3(0.05, 0.12, 0.55)));
	material *blue_200 = new lambertian(new constant_texture(vec3(0.56, 0.79, 0.97)));
	material *light = new diffuse_light(new constant_texture(vec3(45)));
	material *glass = new dielectric(1.3);
	//material *aluminum = new metal(vec3(0.8, 0.85, 0.88), 0.0, "..\\contents\\brdfs\\aluminium.binary");
	//material *aluminimu = new metal(vec3(0.8, 0.85, 0.88), 0.0);
	material *silver = new brdfmaterial("..\\contents\\brdfs\\silver-metallic-paint.binary", vec3(0.8, 0.85, 0.88));
	material *gold_24 = new metal(vec3(0.852, 0.695, 0.449), 0.0);
	material *gold = new metal(vec3(0.945, 0.75, 0.336), 0.0);
	float roughx = 0.02;
	float roughy = 0.3;

	material *beckmann_silver = new beckmann(new constant_texture(vec3(0.8, 0.85, 0.88)), roughx, roughy);
	material *beckmann_white = new beckmann(new constant_texture(vec3(.9, .9, .9)), roughx, roughy);
	material *beckmann_blue = new beckmann(new constant_texture(vec3(0.05, 0.12, 0.55)), roughx, roughy);
	material *beckmann_brow = new beckmann(new constant_texture(vec3(0.426, 0.3, 0.254)), roughx, roughy);
	material *beckmann_gold = new beckmann(new constant_texture(vec3(0.945, 0.75, 0.336)), roughx, roughy);
	//material *beckmann_gold = new beckmann(new constant_texture(vec3(0.852, 0.695, 0.44)), roughx, roughy);

	material* orennayar_white_0 = new orennayar(new constant_texture(vec3(0.7)), 0);
	material* orennayar_blue_200_0 = new orennayar(new constant_texture(vec3(0.56, 0.79, 0.97)), 0);
	material* orennayar_red_200_0 = new orennayar(new constant_texture(vec3(0.93, 0.6, 0.6)), 0);
	material* orennayar_white_10 = new orennayar(new constant_texture(vec3(0.7)), 10);
	material* orennayar_white_20 = new orennayar(new constant_texture(vec3(0.7)), 20);


	//list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, orennayar_blue_200_0));		// right
	//list[i++] = new yz_rect(0, 555, 0, 555, 0, orennayar_red_200_0);							// left
	//list[i++] = new flip_normals(new xz_rect(203, 353, 217, 343, 545, light));	// light
	list[i++] = new flip_normals(new xz_rect(203, 353, 217, 343, 800, light));	// light
	//list[i++] = new sphere(vec3(278, 800, 280), 0.1, light);
	//list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 555, orennayar_white_10));		// top
	list[i++] = new xz_rect(0, 555, 0, 555, 0, orennayar_white_0);							// bottom
	//list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, orennayar_white_10));		// back

	// environment 
	int tx, ty, tn;
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);							// front
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));
	//list[i++] = new env_sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty)));

	model *bunny = new model("..\\contents\\models\\bunny.ply", false, true, orennayar_white_10, vec3(2000, 2000, 2000));
	hitable* b = new translate(new rotate_y(new bvh_node(bunny->genhitablemodel(), bunny->gettrianglecount(), 0, 1), 180), vec3(250, -70, 400));
	//model *dragon = new model("..\\contents\\models\\dragon.ply", false, false, beckmann_gold, vec3(2000, 2000, 2000));
	//model *teapot = new model("..\\contents\\models\\teapot.obj", false, false, beckmann_gold, vec3(80, 80, 80));
	//hitable* b = new translate(new rotate_y(new bvh_node(teapot->genhitablemodel(), teapot->gettrianglecount(), 0, 1), 0), vec3(250, 0, 400));
	list[i++] = b;
	//list[i++] = new constant_medium(b, 0.2, new constant_texture(vec3(0.3, 0.8, 0.2)));

	
	*scene = new hitable_list(list, i);
	

	hitable* light_shape = new flip_normals(new xz_rect(203, 353, 217, 343, 800, 0));
	//hitable* light_shape = new sphere(vec3(278, 800, 280), 0.1, 0);


	//model *bunny_l = new model("..\\contents\\models\\bunny.ply", true, 0, vec3(1000, 1000, 1000));  // bvh_node 的value_pdf 和random是不正确的
	//hitable* b_l = new translate(new rotate_y(bunny->genhitablemodel(), 180), vec3(165, -34, 200));
	hitable** a = new hitable*[7];
	a[0] = light_shape;// new flip_normals(new xz_rect(213, 343, 227, 333, 554, 0));;
	a[1] = new flip_normals(new sphere(lookfrom, 10000, 0));
	//a[1] = new env_sphere(lookfrom, 10000, 0);
	/*a[1] = env_front;
	a[2] = env_back;
	a[3] = env_left;
	a[4] = env_right;
	a[5] = env_top;
	a[6] = env_bottom;*/
	//a[1] = glass_shape;
	//a[1] = b_l;
	*hlist = new hitable_list(a, 1);
}

void teapot_scene(hitable **scene, camera **cam, hitable **hlist, float aspect) {

	vec3 lookfrom(100, 800, -400);
	//vec3 lookfrom(300, 500, -800);
	vec3 lookat(300, 278, 200);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	hitable **list = new hitable*[200];
	
	material *light = new diffuse_light(new constant_texture(vec3(40)));
	float roughx = 0;
	float roughy = 0;

	
	material *beckmann_gold = new beckmann(new constant_texture(vec3(0.945, 0.75, 0.336)), 0.01, 0.05);
	material *beckmann_white = new beckmann(new constant_texture(vec3(1)), roughx, roughy);
	material *beckmann_silver = new beckmann(new constant_texture(vec3(0.8, 0.85, 0.88)), 0.1, 0.1);


	material* lambertian_brow = new lambertian(new constant_texture(vec3(0.426, 0.3, 0.254)));
	material* orennayar_brow = new orennayar(new constant_texture(vec3(0.426, 0.3, 0.254)), 0);
	material* orennayar_white = new orennayar(new constant_texture(vec3(1)), 10);
	
	material* glass = new dielectric(1.3);


	list[i++] = new flip_normals(new xz_rect(3, 153, 217, 343, 800, light));	// light
	//list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 1500, light));	// light
	list[i++] = new xz_rect(0, 555, 0, 555, 0, lambertian_brow);							// bottom

	// environment 
	int tx, ty, tn;
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);							// front
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));

	/*model *teapot = new model("..\\contents\\models\\teapot.obj", false, false, beckmann_gold, vec3(80, 80, 80));
	hitable* b = new translate(new rotate_y(new bvh_node(teapot->genhitablemodel(), teapot->gettrianglecount(), 0, 1), 0), vec3(250, 0, 400));
	list[i++] = b;*/
	teapot* newTeapot = new teapot(40, new metal(vec3(0.9), 0.0));
	hitable** teapotList = newTeapot->createPloyTeapot();
	list[i++] = new translate(new rotate_x(new bvh_node(teapotList, newTeapot->getTriangleCount(), 0, 1), 90), vec3(200, 0, 250));

	teapot* newTeapot2 = new teapot(40, beckmann_gold/*new metal(vec3(0.9), 0.0)*/);
	hitable** teapotList2 = newTeapot2->createPloyTeapot();
	list[i++] = new translate(new rotate_x(new bvh_node(teapotList2, newTeapot2->getTriangleCount(), 0, 1), 90), vec3(360, 0, 150));

	model *bunny = new model("..\\contents\\models\\bunny.ply", false, true, orennayar_white, vec3(2000, 2000, 2000));
	hitable* b = new translate(new rotate_y(new bvh_node(bunny->genhitablemodel(), bunny->gettrianglecount(), 0, 1), 180), vec3(180, -70, 450));
	list[i++] = b;
	
	hitable* c = new sphere(vec3(280, 30, 70), 30, new dielectric(1.5));
	list[i++] = c;
	list[i++] = new constant_medium(c, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));

	model *dragon = new model("..\\contents\\models\\dragon.ply", false, true, beckmann_silver, vec3(500, 500, 500));
	hitable* d = new translate(new rotate_y(new bvh_node(dragon->genhitablemodel(), dragon->gettrianglecount(), 0, 1), 180), vec3(140, -20, 120));
	list[i++] = d;

	*scene = new hitable_list(list, i);


	hitable* light_shape = new flip_normals(new xz_rect(3, 153, 217, 343, 800, 0));


	hitable** a = new hitable*[1];
	a[0] = light_shape;// new flip_normals(new xz_rect(213, 343, 227, 333, 554, 0));;
	*hlist = new hitable_list(a, 1);
}

void ball_scenes(hitable **scene, camera **cam, hitable **hlist, float aspect) {

	vec3 lookfrom(300, 600, -100);
	vec3 lookat(300, 20, 250);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	hitable **list = new hitable*[200];
	
	material *light = new diffuse_light(new constant_texture(vec3(20)));

	material* orennayar_white_10 = new orennayar(new constant_texture(vec3(0.7)), 10);
	material* orennayar_brow = new orennayar(new constant_texture(vec3(0.426, 0.3, 0.254)), 0);

	list[i++] = new flip_normals(new xz_rect(203, 353, 217, 343, 800, light));	// light
	list[i++] = new xz_rect(-100, 655, -100, 655, 0, orennayar_brow);							// bottom

	// environment 
	int tx, ty, tn;
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);							// front
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));
	
	material **mats = new material*[121];
	for (size_t j = 0; j < 121; j++)
	{

		mats[j] = new beckmann(new constant_texture(vec3(1)), (j % 11) / 11.0f, (j / 11) / 11.0f);
		//mats[j] = new beckmann(new constant_texture(vec3(0.945, 0.75, 0.336)), (j % 11) / 11.0f, (j / 11) / 11.0f);
		list[i++] = new sphere(vec3(550 - (j % 11) * 50, 20, 450 - 50 * (j / 11)), 20, mats[j]);
	}

	
	*scene = new hitable_list(list, i);


	hitable* light_shape = new flip_normals(new xz_rect(203, 353, 217, 343, 800, 0));



	hitable** a = new hitable*[1];
	a[0] = light_shape;
	
	*hlist = new hitable_list(a, 1);
}

void ball_orennayar_scenes(hitable **scene, camera **cam, hitable **hlist, float aspect) {

	vec3 lookfrom(300, 800, -100);
	vec3 lookat(300, 20, 450);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	hitable **list = new hitable*[21];

	material *light = new diffuse_light(new constant_texture(vec3(20)));

	material* orennayar_white_10 = new orennayar(new constant_texture(vec3(0.7)), 10);
	material* orennayar_brow = new orennayar(new constant_texture(vec3(0.426, 0.3, 0.254)), 0);

	list[i++] = new flip_normals(new xz_rect(203, 353, 217, 343, 800, light));	// light
	list[i++] = new xz_rect(-100, 655, -100, 655, 0, orennayar_brow);							// bottom

	// environment 
	int tx, ty, tn;
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);							// front
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));

	material **mats = new material*[21];
	for (size_t j = 0; j <= 20; j++)
	{

		mats[j] = new orennayar(new constant_texture(vec3(1)), j);
		//mats[j] = new beckmann(new constant_texture(vec3(0.945, 0.75, 0.336)), (j % 11) / 11.0f, (j / 11) / 11.0f);
		list[i++] = new sphere(vec3(550 - (j % 7) * 70, 30, 450 - 70 * (j / 7)), 30, mats[j]);
	}


	*scene = new hitable_list(list, i);


	hitable* light_shape = new flip_normals(new xz_rect(203, 353, 217, 343, 800, 0));



	hitable** a = new hitable*[1];
	a[0] = light_shape;

	*hlist = new hitable_list(a, 1);
}

void final(hitable **scene, camera **cam, hitable **hlist, float aspect) {
	int nb = 20;
	hitable **list = new hitable*[30];
	hitable **boxlist = new hitable*[10000];
	hitable **boxlist2 = new hitable*[10000];
	material *white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
	material *ground = new lambertian(new constant_texture(vec3(0.48, 0.83, 0.53)));
	int b = 0;
	for (int i = 0; i < nb; i++) {
		for (int j = 0; j < nb; j++) {
			float w = 100;
			float x0 = -1000 + i * w;
			float z0 = -1000 + j * w;
			float y0 = 0;
			float x1 = x0 + w;
			float y1 = 100 * (drand48() + 0.01);
			float z1 = z0 + w;
			boxlist[b++] = new box(vec3(x0, y0, z0), vec3(x1, y1, z1), ground);
		}
	}
	int l = 0;
	list[l++] = new bvh_node(boxlist, b, 0, 1);
	material *light = new diffuse_light(new constant_texture(vec3(7, 7, 7)));
	list[l++] = new flip_normals(new xz_rect(123, 423, 147, 412, 554, light));
	vec3 center(400, 400, 200);
	list[l++] = new moving_sphere(center, center + vec3(30, 0, 0), 0, 1, 50, new lambertian(new constant_texture(vec3(0.7, 0.3, 0.1))));
	list[l++] = new sphere(vec3(260, 150, 45), 50, new dielectric(1.5));
	list[l++] = new sphere(vec3(0, 150, 145), 50, new metal(vec3(0.8, 0.8, 0.9), 1.0));// , 10.0));
	hitable *boundary = new sphere(vec3(360, 150, 145), 70, new dielectric(1.5));
	list[l++] = boundary;
	list[l++] = new constant_medium(boundary, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));
	boundary = new sphere(vec3(0, 0, 0), 5000, new dielectric(1.5));
	list[l++] = new constant_medium(boundary, 0.0001, new constant_texture(vec3(1.0, 1.0, 1.0)));
	int nx, ny, nn;
	unsigned char *tex_data = stbi_load("..\\contents\\textures\\earthmap.jpg", &nx, &ny, &nn, 0);
	material *emat = new lambertian(new image_texture(tex_data, nx, ny));
	list[l++] = new sphere(vec3(400, 200, 400), 100, emat);
	texture *pertext = new noise_texture(0.1);
	list[l++] = new sphere(vec3(220, 280, 300), 80, new lambertian(pertext));
	int ns = 1000;
	for (int j = 0; j < ns; j++) {
		boxlist2[j] = new sphere(vec3(165 * drand48(), 165 * drand48(), 165 * drand48()), 10, white);
	}
	list[l++] = new translate(new rotate_y(new bvh_node(boxlist2, ns, 0.0, 1.0), 15), vec3(-100, 270, 395));

	*scene = new hitable_list(list, l);
	vec3 lookfrom(478, 278, -600);
	vec3 lookat(278, 278, 0);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	hitable* light_shape = new flip_normals(new xz_rect(123, 423, 147, 412, 554, 0));
	hitable** a = new hitable*[1];
	a[0] = light_shape;

	*hlist = new hitable_list(a, 1);
}

void jadebunny_scene(hitable **scene, camera ** cam, hitable **hlist, float aspect)
{
	vec3 lookfrom(300, 500, -800);
	vec3 lookat(300, 278, 200);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	brdf* aluminum_brdf = new brdf();
	hitable **list = new hitable*[200];
	material *light = new diffuse_light(new constant_texture(vec3(45)));
	material *glass = new dielectric(1.2);


	material* orennayar_white_0 = new orennayar(new constant_texture(vec3(0.7)), 0);
	material* orennayar_blue = new orennayar(new constant_texture(vec3(0.2, 0.4, 0.9)), 0);


	list[i++] = new flip_normals(new xz_rect(203, 353, 17, 543, 800, light));	// light
	list[i++] = new xz_rect(0, 555, 0, 555, 0, orennayar_white_0);							// bottom

	// environment 
	int tx, ty, tn;
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);							// front
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));

	model *bunny = new model("..\\contents\\models\\bunny.ply", false, false, glass, vec3(2000, 2000, 2000));
	hitable* b = new translate(new rotate_y(new bvh_node(bunny->genhitablemodel(), bunny->gettrianglecount(), 0, 1), 180), vec3(250, -70, 400));
	list[i++] = b;
	//vec3(0.2, 0.4, 0.9) 0.3, 0.8, 0.2
	//list[i++] = new constant_medium(b, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));
	model *bunny1 = new model("..\\contents\\models\\bunny.ply", false, true, orennayar_blue, vec3(1990, 1990, 1990));
	//hitable* c = new translate(new rotate_y(new bvh_node(bunny1->genhitablemodel(), bunny1->gettrianglecount(), 0, 1), 180), vec3(250, -70, 400));
	list[i++] = new translate(new rotate_y(new bvh_node(bunny1->genhitablemodel(), bunny1->gettrianglecount(), 0, 1), 180), vec3(250, -70, 400));
	//list[i++] = new constant_medium(c, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));


	*scene = new hitable_list(list, i);


	hitable* light_shape = new flip_normals(new xz_rect(203, 353, 217, 343, 800, 0));


	hitable** a = new hitable*[7];
	a[0] = light_shape;// new flip_normals(new xz_rect(213, 343, 227, 333, 554, 0));;
	*hlist = new hitable_list(a, 1);
}

void soldier_scene(hitable **scene, camera ** cam, hitable **hlist, float aspect)
{
	vec3 lookfrom(300, 500, -800);
	vec3 lookat(300, 278, 200);
	float dist_to_focus = 1000.0;
	float aperture = 10;
	float vfov = 40.0;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	brdf* aluminum_brdf = new brdf();
	hitable **list = new hitable*[200];
	material *light = new diffuse_light(new constant_texture(vec3(45)));
	material *light_1 = new diffuse_light(new constant_texture(vec3(35)));


	material* orennayar_white_0 = new orennayar(new constant_texture(vec3(0.7)), 0);
	material* orennayar_blue = new orennayar(new constant_texture(vec3(0.2, 0.4, 0.9)), 0);

	float roughx = 0.9;
	float roughy = 0.85;
	material *beckmann_silver = new beckmann(new constant_texture(vec3(0.8, 0.85, 0.88)), roughx, roughy);
	material *beckmann_white = new beckmann(new constant_texture(vec3(.9, .9, .9)), 0, 0);
	material *beckmann_blue = new beckmann(new constant_texture(vec3(0.05, 0.12, 0.55)), roughx, roughy);
	material *beckmann_brow = new beckmann(new constant_texture(vec3(0.426, 0.3, 0.254)), roughx, roughy);
	material *beckmann_gold = new beckmann(new constant_texture(vec3(0.945, 0.75, 0.336)), roughx, roughy);

	material *beckmann_bottom = new beckmann(new constant_texture(vec3(1)), 0.9, 0.9);
	int tx, ty, tn;
	unsigned char* tex_floor = stbi_load("..\\contents\\textures\\TexturesCom_Wood_Wenge_1K_albedo.png", &tx, &ty, &tn, 0);							// 

	//
	material *glass = new dielectric(1.4);

	material *imaged_bottom = new orennayar(new image_texture(tex_floor, tx, ty), 0.5);

	//list[i++] = new flip_normals(new xz_rect(203, 353, 17, 543, 800, light));	// light
	//list[i++] = new xz_rect(0, 555, 0, 555, 0, orennayar_white_0);							// bottom
	list[i++] = new flip_normals(new xz_rect(203, 353, 17, 167, 800, light_1));	// light
	//list[i++] = new xz_rect(0, 600, 0, 600, 0, imaged_bottom);							// bottom
	//list[i++] = new xz_rect(0, 600, 0, 600, 1, glass);							// bottom2
	list[i++] = new box(vec3(0, -0.1, 0), vec3(600, 0.1, 600), imaged_bottom);

	list[i++] = new box(vec3(0, -1, 0), vec3(600, 1, 600), glass);
	// environment 
	//unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);							// front
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky4.jpg", &tx, &ty, &tn, 0);							// front
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));

	int nx, ny, nn;
	//unsigned char* tex_soldier = stbi_load("..\\contents\\textures\\NPC_ChaoJiBing_A.png", &nx, &ny, &nn, 0);
	unsigned char* tex_soldier = stbi_load("..\\contents\\textures\\NPC_YuanChengBing_A.png", &nx, &ny, &nn, 0);
	//unsigned char* tex_soldier = stbi_load("..\\contents\\textures\\hero08_skin1_D.png", &nx, &ny, &nn, 0);
	material *beckmann_tex = new beckmann(new image_texture(tex_soldier, nx, ny), roughx, roughy);
	//model *soilder = new model("..\\contents\\models\\Soilder.FBX", false, true, beckmann_tex, vec3(2.5f, 2.5f, 2.5f));
	model *soilder = new model("..\\contents\\models\\Soilder.FBX", false, true, beckmann_tex, vec3(8, 8, 8));
	//hitable* b = new translate(new rotate_x(new rotate_y(new bvh_node(soilder->genhitablemodel(), soilder->gettrianglecount(), 0, 1), 180), -90), vec3(250, 0, 300));
	hitable* b = new translate(new rotate_y(new bvh_node(soilder->genhitablemodel(), soilder->gettrianglecount(), 0, 1), 180), vec3(250, 0, 300));
	list[i++] = b;

	*scene = new hitable_list(list, i);


	hitable* light_shape = new flip_normals(new xz_rect(203, 353, 17, 167, 800, 0));
	hitable* mirror_shape = new rotate_y(new xy_rect(0, 400, 0, 800, 800, 0), 30);


	hitable** a = new hitable*[7];
	a[0] = light_shape;
	a[1] = mirror_shape;
	*hlist = new hitable_list(a, 1);

}

void flatnormal_bunny(hitable **scene, camera **cam, hitable **hlist, float aspect)
{
	vec3 lookfrom(300, 500, -800);
	vec3 lookat(300, 278, 200);
	float dist_to_focus = 10.0;
	float aperture = 0;
	float vfov = 40;
	*cam = new camera(lookfrom, lookat, vec3(0, 1, 0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);

	int i = 0;
	hitable **list = new hitable*[200];
	material *light = new diffuse_light(new constant_texture(vec3(45)));

	material* orennayar_white = new orennayar(new constant_texture(vec3(0.7)), 0.1f);
	material* beckmann_gold = new beckmann(new constant_texture(vec3(0.945, 0.75, 0.336)), 0.85f, 0.85f);

	list[i++] = new flip_normals(new xz_rect(203, 353, 17, 167, 800, light));
	list[i++] = new xz_rect(0, 600, 0, 600, 0, orennayar_white);

	// environment
	int tx, ty, tn;
	unsigned char* tex_data = stbi_load("..\\contents\\environment_map\\sky_2.png", &tx, &ty, &tn, 0);
	list[i++] = new flip_normals(new sphere(lookfrom, 10000, new diffuse_light(new image_texture(tex_data, tx, ty))));

	model *bunny = new model("..\\contents\\models\\bunny.ply", false, false, beckmann_gold, vec3(2000));
	hitable* b = new translate(new rotate_y(new bvh_node(bunny->genhitablemodel(), bunny->gettrianglecount(), 0, 1), 180), vec3(250, -70, 400));

	*scene = new hitable_list(list, i);

	hitable* light_shape = new flip_normals(new xz_rect(203, 353, 17, 167, 800, 0));
	hitable** a = new hitable*[7];
	a[0] = light_shape;
}

hitable* final1()
{
	int nb = 20;
	hitable **list = new hitable*[10000];
	hitable **boxlist = new hitable*[10000];
	hitable **boxlist2 = new hitable*[10000];
	material *white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
	int l = 0;
	material *light = new diffuse_light(new constant_texture(vec3(7, 7, 7)));
	list[l++] = new xz_rect(123, 423, 147, 412, 554, light);
	vec3 center(400, 400, 200);
	int ns = 1000;
	for (int j = 0; j < ns; j++) {
		//list[l++] = new translate(new sphere(vec3(165 * drand48(), 165 * drand48(), 165 * drand48()), 10, white), vec3(-100, 270, 395));
		boxlist2[j] = new sphere(vec3(165 * drand48(), 165 * drand48(), 165 * drand48()), 10, white);
	}
	list[l++] = new translate(new rotate_y(new bvh_node(boxlist2, ns, 0.0, 1.0), 15), vec3(-100, 270, 395));
	return new hitable_list(list, l);
}


bool isfinished = false;

int getpixel()
{
	return finishedPixel++;
}

double **sobol_points(unsigned N, unsigned D, const char *dir_file)
{
	ifstream infile(dir_file, ios::in);
	if (!infile) {
		cout << "Input file containing direction numbers cannot be found!\n";
		exit(1);
	}
	char buffer[1000];
	infile.getline(buffer, 1000, '\n');

	// L = max number of bits needed 
	unsigned L = (unsigned)ceil(log((double)N) / log(2.0));

	// C[i] = index from the right of the first zero bit of i
	unsigned *C = new unsigned[N];
	C[0] = 1;
	for (unsigned i = 1; i <= N - 1; i++) {
		C[i] = 1;
		unsigned value = i;
		while (value & 1) {
			value >>= 1;
			C[i]++;
		}
	}

	// POINTS[i][j] = the jth component of the ith point
	//                with i indexed from 0 to N-1 and j indexed from 0 to D-1
	double **POINTS = new double *[N];
	for (unsigned i = 0; i <= N - 1; i++) POINTS[i] = new double[D];
	for (unsigned j = 0; j <= D - 1; j++) POINTS[0][j] = 0;

	// ----- Compute the first dimension -----

	// Compute direction numbers V[1] to V[L], scaled by pow(2,32)
	unsigned *V = new unsigned[L + 1];
	for (unsigned i = 1; i <= L; i++) V[i] = 1 << (32 - i); // all m's = 1

	// Evalulate X[0] to X[N-1], scaled by pow(2,32)
	unsigned *X = new unsigned[N];
	X[0] = 0;
	for (unsigned i = 1; i <= N - 1; i++) {
		X[i] = X[i - 1] ^ V[C[i - 1]];
		POINTS[i][0] = (double)X[i] / pow(2.0, 32); // *** the actual points
		//        ^ 0 for first dimension
	}

	// Clean up
	delete[] V;
	delete[] X;


	// ----- Compute the remaining dimensions -----
	for (unsigned j = 1; j <= D - 1; j++) {

		// Read in parameters from file 
		unsigned d, s;
		unsigned a;
		infile >> d >> s >> a;
		unsigned *m = new unsigned[s + 1];
		for (unsigned i = 1; i <= s; i++) infile >> m[i];
		// Compute direction numbers V[1] to V[L], scaled by pow(2,32)
		unsigned *V = new unsigned[L + 1];
		if (L <= s) {
			for (unsigned i = 1; i <= L; i++) V[i] = m[i] << (32 - i);
		}
		else {
			for (unsigned i = 1; i <= s; i++) V[i] = m[i] << (32 - i);
			for (unsigned i = s + 1; i <= L; i++) {
				V[i] = V[i - s] ^ (V[i - s] >> s);
				for (unsigned k = 1; k <= s - 1; k++)
					V[i] ^= (((a >> (s - 1 - k)) & 1) * V[i - k]);
			}
		}

		// Evalulate X[0] to X[N-1], scaled by pow(2,32)
		unsigned *X = new unsigned[N];
		X[0] = 0;
		for (unsigned i = 1; i <= N - 1; i++) {
			X[i] = X[i - 1] ^ V[C[i - 1]];
			POINTS[i][j] = (double)X[i] / pow(2.0, 32); // *** the actual points
			//        ^ j for dimension (j+1)
		}

		// Clean up
		delete[] m;
		delete[] V;
		delete[] X;
	}
	delete[] C;

	return POINTS;
}


void renderthread(hitable* world, hitable* hlist, camera *cam, double** samplepoints)
{
	while (finishedPixel < nx*ny)
	{
		g_lock.lock();
		int curpixel = finishedPixel++;
		if ((finishedPixel - 1) % ny == 0)
		{
			std::cout << "\r............" << (float)curpixel / (float)(nx*ny) * 100 << "%" << std::flush;
		}
		g_lock.unlock();
		vec3 col(0, 0, 0);
		int i = curpixel % ny;
		int j = ny - 1 - (curpixel / nx);
		int goodsample_count = 0;
		for (int s = 0; s < ns; s++)
		{
			/*float u = float(i + drand48()) / float(nx);
			float v = float(j + drand48()) / float(ny);*/
			float u = float(samplepoints[s][0] + i) / float(nx);
			float v = float(samplepoints[s][1] + j) / float(ny);
			ray r = cam->get_ray(u, v);
			//vec3 p = r.point_at_parameter(2.0);
			int* depth = (int*)malloc(sizeof(int));
			*depth = 0;
			//cout << *depth << endl;
			col += de_nan(color(r, world, hlist, depth));
			free(depth);
			if (*depth <= maxDepth)
			{
				goodsample_count++;
			}
		}
		col /= float(ns);

		col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
		/*if (col.length() > 1)
			col.make_unit_vector();*/
		int ir = int(255.99 * col[0]);
		int ig = int(255.99 * col[1]);
		int ib = int(255.99 * col[2]);

		ir = ir > 255 ? 255 : ir;
		ig = ig > 255 ? 255 : ig;
		ib = ib > 255 ? 255 : ib;
		ir = ir < 0 ? 0 : ir;
		ig = ig < 0 ? 0 : ig;
		ib = ib < 0 ? 0 : ib;
		int index = curpixel;
		colors[index] = new int[3];
		colors[index][0] = ir;
		colors[index][1] = ig;
		colors[index][2] = ib;
	}
	g_lock.lock();
	if (!isfinished)
	{
		isfinished = true;
		for (int i = 0; i < nx*ny; i++)
		{
			outfile << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << "\n";
		}
	}
	g_lock.unlock();
}


int main()
{
#ifdef RaysBackgroundY
	//ofstream outfile("..\\results\\final.ppm", ios_base::out);
	outfile << "P3\n" << nx << " " << ny << "\n255\n";

	std::cout << "P3\n" << nx << " " << ny << "\n255\n";

	hitable *world;
	camera *cam;
	hitable *hlist;

	switch (sceneid)
	{
	case 0:
		cornell_box(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 1:
		teapot_scene(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 2:
		ball_scenes(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 3:
		ball_orennayar_scenes(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 4:
		jadebunny_scene(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 5:
		final(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 6:
		soldier_scene(&world, &cam, &hlist, float(nx) / float(ny));
		break;
	case 7:
		flatnormal_bunny(&world, &cam, &hlist, float(nx) / float(ny));
	default:
		break;
	}
	

	int total = nx * ny;

	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
	int startmscout = ms.count();
	colors = new int*[ny*nx];
	
	double** samplepoints = sobol_points(ns, 2, "..\\contents\\sobol\\new-joe-kuo-6.21201");

	thread t[thread_count];
	for (int i = 0; i < thread_count; i++)
	{
		t[i] = thread(renderthread, world, hlist, cam, samplepoints);
	}
	curfinishedcount = 0;
	for (int i = 0; i < thread_count; i++)
	{
		t[i].join();
	}
	
	std::cout << endl;
	ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
	int endmscount = ms.count();
	std::cout << endmscount - startmscout<< "ms" << endl;

#endif // RaysBackgroundY

	_CrtDumpMemoryLeaks();
	return 0;
}

