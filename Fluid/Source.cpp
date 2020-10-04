#define SFML_STATIC
#include <SFML/Graphics.hpp>
#include <glm/glm.hpp>
#include <iostream>

const float scale = 5.0f;

const int WIDTH = 800/scale;
const int HEIGHT = 400/scale;
const float c = 1;

const float step = 5;
const int numVertex = WIDTH / 5;

std::vector<glm::vec3> dyes;
std::vector<glm::vec3> newDyes;
std::vector<glm::vec2> oldVelocities;
std::vector<glm::vec2> velocities;
std::vector<glm::vec2> newVelocities;
std::vector<float> pressure;
std::vector<bool> isSolid;

float dt = 0.01f;


float randf() {
	return ((float)rand() / RAND_MAX) * 2 - 1;
}

float randf01() {
	return ((float)rand() / RAND_MAX);
}

void createDye(int sizeX, int sizeY) {
	dyes.reserve(sizeX * sizeY);
	newDyes.reserve(sizeX * sizeY);
	for (int i = 0; i < (sizeX * sizeY); i++) {
		int x = i % WIDTH;
		int y = i / WIDTH;

		if(x>WIDTH/2 && x <3*WIDTH/4 && y > HEIGHT/2 && y < 3* HEIGHT/4)
			dyes.emplace_back(randf01(), 0.1f, 0.4f);
		else
			dyes.emplace_back(0, 0.0f, 0.0f);
		newDyes.push_back(dyes[i]);
	}
}



void createVelocities(int sizeX, int sizeY) {
	velocities.reserve(sizeX * sizeY);
	newVelocities.reserve(sizeX * sizeY);
	for (int y = -sizeY / 2; y < (sizeY / 2); y++) {
		for (int x = -sizeX / 2; x < (sizeX / 2); x++) {
			newVelocities.emplace_back(0);
			if (x == -sizeX / 2 || x == sizeX / 2 - 1 || y == sizeY / 2 - 1 || y == -sizeX / 2)
				velocities.emplace_back(0, 0);
			else
				//velocities.emplace_back((float)(y) / (float)HEIGHT, -(float)(x) / (float)WIDTH);
				velocities.emplace_back(randf()/100.0f, randf()/100.0f);
				//velocities.emplace_back(randf() + (float)(y) / (float)HEIGHT * 4, randf() - (float)(x) / (float)WIDTH * 4);
				//velocities.emplace_back(-(float)(x) / (float)WIDTH, 0);
				//velocities.emplace_back(0,-(float)(y) / (float)HEIGHT);
				//velocities.emplace_back(x < 0 ? 1 : -1, 0);

		}
	}
}

int at(int x, int y) {
	return y * WIDTH + x;
}

float divergence(int x, int y, const std::vector<glm::vec2>& vel) {
	float div = 0;
	div += vel[at(x - 1, y)].x - vel[at(x + 1, y)].x;
	div += vel[at(x, y - 1)].y - vel[at(x, y + 1)].y;
	return div;
}

glm::vec2 divergence2(int x, int y, const std::vector<glm::vec2>& vel) {
	glm::vec2 div = { 0,0 };
	for (int i = -1; i < 2; i++) {
		for (int j = -1; j < 2; j++) {
			div += vel[at(x + i, y + j)] - vel[at(x, y)];
		}
	}
	return div;
}


float divergence21(int x, int y, const std::vector<glm::vec2>& vel) {
	/*
	(dot(f[x-1,y-1] + f[x+1,y+1], {1,1}) +
dot(f[x - 1, y + 1] + f[x + 1, y - 1], { 1,-1 })* { 1, -1 } +
(f[x - 1, y] + f[x + 1, y] - f[x, y - 1] - f[x, y + 1])* {
	2, -2
} +
f[x, y] * -4) * 1 / 8
*/

	glm::vec2 div = { 0,0 };
	auto diag = glm::dot(vel[at(x - 1, y - 1)] + vel[at(x + 1, y + 1)], glm::vec2(1, 1))
		+ glm::dot(vel[at(x - 1, y + 1)] + vel[at(x + 1, y - 1)], { 1,-1 });
	auto horizontal = (vel[at(x - 1, y)].x - vel[at(x + 1, y)].x) * 2.0f;
	auto vertical = (vel[at(x, y - 1)].y - vel[at(x - 1, y)].y) * 2.0f;
	auto mid = vel[at(x, y)] * 4.0f;


	return diag + horizontal + vertical + mid.x + mid.y;
}

void lin_solve(std::vector<float>& x, std::vector<float>& x0, float a, float c, int iter)
{
	float cRecip = 1.0 / c;
	for (int k = 0; k < iter; k++) {
		for (int j = 1; j < HEIGHT - 1; j++) {
			for (int i = 1; i < WIDTH - 1; i++) {
				x[at(i, j)] =
					(x0[at(i, j)]
						+ a * (x[at(i + 1, j)]
							+ x[at(i - 1, j)]
							+ x[at(i, j + 1)]
							+ x[at(i, j - 1)]
							)) * cRecip;
			}
		}
	}
}

std::vector<float> calcDivergence() {
	std::vector<float> div;
	div.reserve(WIDTH * HEIGHT);
	for (int j = 0; j < HEIGHT; j++) {
		for (int i = 0; i < WIDTH; i++) {
			if (i == 0 || j == 0 || i == WIDTH - 1 || j == HEIGHT - 1) {
				div.emplace_back(0);
			}
			else {
				div.emplace_back(
					-0.5f * (
						velocities[at(i + 1, j)].x
						- velocities[at(i - 1, j)].x
						+ velocities[at(i, j + 1)].y
						- velocities[at(i, j - 1)].y
						) );
			}
		}
	}
	return div;
}

std::vector<float> createPressure() {
	std::vector<float> pressure;
	pressure.reserve(WIDTH * HEIGHT);
	for (int j = 0; j < HEIGHT; j++) {
		for (int i = 0; i < WIDTH; i++) {
			pressure.emplace_back(0);
		}
	}
	return pressure;
}



void divergentVelocity() {

	auto divergences = calcDivergence();
	auto pressures = createPressure();
	int iter = 10;
	lin_solve(pressures, divergences, 1, 6, 10);

	for (int y = 1; y < HEIGHT - 1; y++) {
		for (int x = 1; x < WIDTH - 1; x++) {

			float dj = 0.001f;

			/*float ddiv_x_dt = divergence(x-1, y, velocities) - divergence(x+1, y, velocities);
			float ddiv_y_dt = divergence(x, y-1, velocities) - divergence(x, y+1, velocities);


			velocities[index].x += ddiv_x_dt*dj;
			velocities[index].y += ddiv_y_dt*dj;*/

			/*auto div = divergence2(x, y, velocities);
			velocities[at(x,y)] += div * dj;
			auto div1 = divergence2(x, y, velocities);
			int d = 0;*/

			/*auto div = divergence21(x, y, velocities);
			velocities[at(x, y)] -= velocities[at(x, y)]* div * dj;
			auto div1 = divergence21(x, y, velocities);
			int d = 0;*/


			velocities[at(x, y)].x -= 0.5f * (pressures[at(x + 1, y)] - pressures[at(x - 1, y)]);
			velocities[at(x, y)].y -= 0.5f * (pressures[at(x, y + 1)] - pressures[at(x, y - 1)]);
		}
	}
}

//void createBoundary(int sizeX, int sizeY) {
//	bounda.reserve(sizeX * sizeY);
//	for (int y = -sizeY / 2; y < (sizeY / 2); y++) {
//		for (int x = -sizeX / 2; x < (sizeX / 2); x++) {
//			//velocities.emplace_back(0, 0);
//			velocities.emplace_back((float)(y) / (float)HEIGHT, -(float)(x) / (float)WIDTH);
//		}
//	}
//}


int roundi(float x) {
	x = x + 0.5 - (x < 0);
	return (int)x;
}

void moveTo(glm::vec2 dest, glm::vec3 value, int i) {
	int destx = roundi(dest.x);
	int desty = roundi(dest.y);
	if (desty > 0 && desty < HEIGHT && destx > 0 && destx < WIDTH) {
		int destIndex = desty * WIDTH + destx;
		newDyes[destIndex] += value;
	}
	else {
		newDyes[i] += value;
	}
}



void convect() {
	for (int i = 0; i < velocities.size(); i++) {
		auto dye = dyes[i];
		auto vel = velocities[i];
		int x = i % WIDTH;
		//11  21

		int y = i / WIDTH;
		auto pos = glm::vec2(x, y);
		auto dest = pos + vel * dt;

		//21  22
		//dest -= glm::vec2(0.5f, 0.5f);

		auto desti = glm::vec2(roundi(dest.x), roundi(dest.y));
		float distances[3][3];
		float total = 0;
		for (int k = -1; k <= 1; k++) {
			for (int j = -1; j <= 1; j++) {
				auto destxy = desti + glm::vec2(k, j);
				distances[k + 1][j + 1] = 0;
				if (destxy.y > 0 && destxy.y < HEIGHT && destxy.x > 0 && destxy.x < WIDTH) {
					float dif_x = abs(destxy.x - dest.x);
					float dif_y = abs(destxy.y - dest.y);
					if (dif_x < 1 && dif_y < 1) {
						float dis = glm::distance(destxy, dest);
						distances[k + 1][j + 1] = 1 - dis;
						total += distances[k + 1][j + 1];
					}
				}
			}
		}

		for (int k = -1; k <= 1; k++) {
			for (int j = -1; j <= 1; j++) {
				auto destxy = desti + glm::vec2(k, j);
				if (destxy.y > 0 && destxy.y < HEIGHT && destxy.x > 0 && destxy.x < WIDTH) {
					moveTo(destxy, dye * distances[k + 1][j + 1] / total, i);
				}
			}
		}



		//std::cout << "pos " << x << "," << y << std::endl;
		//std::cout << "p11 " << p11 << " " << roundi(dest11.x) << "," << roundi(dest11.y) << std::endl;
		//std::cout << "p12 " << p12 << " " << roundi(dest12.x) << "," << roundi(dest12.y) << std::endl;
		//std::cout << "p21 " << p21 << " " << roundi(dest21.x) << "," << roundi(dest21.y) << std::endl;
		//std::cout << "p22 " << p22 << " " << roundi(dest22.x) << "," << roundi(dest22.y) << std::endl;

		////at dest11
		//moveTo(dest11, dye*p11, i);
		//moveTo(dest12, dye*p12, i);
		//moveTo(dest21, dye*p21, i);
		//moveTo(dest22, dye*p22, i);
	}
}

template <typename T>
T linearInterpolation(float x1, T f_x1, float x2, T f_x2, float x)
{
	T result = f_x2 * ((x - x1) / (x2 - x1)) + f_x1 * ((x2 - x) / (x2 - x1));
	return result;
}

void convectback() {
	for (int i = 0; i < velocities.size(); i++) {
		auto vel = velocities[i];
		int x = i % WIDTH;
		int y = i / WIDTH;

		auto pos = glm::vec2(x, y);
		auto dest = pos - vel * dt;


		float x1 = floor(dest.x);
		float x2 = x1 + 1;
		float y1 = floor(dest.y);
		float y2 = y1 + 1;

		int xi = (int)dest.x;
		int yi = (int)dest.y;

		if (yi >= 0 && yi < HEIGHT - 1 && xi >= 0 && xi < WIDTH - 1) {
			auto R1 = linearInterpolation(x1, dyes[at(xi, yi)], x2, dyes[at(xi + 1, yi)], dest.x);
			auto R2 = linearInterpolation(x1, dyes[at(xi, yi + 1)], x2, dyes[at(xi + 1, yi + 1)], dest.x);
			auto  P = linearInterpolation(y1, R1, y2, R2, dest.y);


			newDyes[at(x, y)] = P;


			auto V1 = linearInterpolation(x1, velocities[at(xi, yi)], x2, velocities[at(xi + 1, yi)], dest.x);
			auto V2 = linearInterpolation(x1, velocities[at(xi, yi + 1)], x2, velocities[at(xi + 1, yi + 1)], dest.x);
			auto V3 = linearInterpolation(y1, V1, y2, V2, dest.y);


			newVelocities[at(x, y)] = V3;
		}

	}
}




int main()
{
	glm::vec2 dest = glm::vec2(1, 1);
	dest -= glm::vec2(0.5f, 0.5f);

	auto dest11 = glm::vec2(floor(dest.x), floor(dest.y));
	auto dest12 = glm::vec2(floor(dest.x), ceil(dest.y));;
	auto dest21 = glm::vec2(ceil(dest.x), floor(dest.y));;
	auto dest22 = glm::vec2(ceil(dest.x), ceil(dest.y));;

	auto p11 = (dest.x - dest11.x) * (dest.y - dest11.y);
	auto p12 = (dest.x - dest12.x) * -(dest.y - dest12.y);
	auto p21 = -(dest.x - dest21.x) * (dest.y - dest21.y);
	auto p22 = -(dest.x - dest22.x) * -(dest.y - dest22.y);

	std::cout << "p11 " << p11 << " " << roundi(dest11.x) << "," << roundi(dest11.y) << std::endl;
	std::cout << "p12 " << p12 << " " << roundi(dest12.x) << "," << roundi(dest12.y) << std::endl;
	std::cout << "p21 " << p21 << " " << roundi(dest21.x) << "," << roundi(dest21.y) << std::endl;
	std::cout << "p22 " << p22 << " " << roundi(dest22.x) << "," << roundi(dest22.y) << std::endl;







	sf::RenderWindow window(sf::VideoMode(800, 400), "SFML works!");



	createDye(WIDTH, HEIGHT);
	createVelocities(WIDTH, HEIGHT);


	sf::Sprite fluidSprite;
	fluidSprite.scale(sf::Vector2f(scale, scale));
	sf::Texture fluidTexture;
	
	sf::Image image;
	fluidTexture.create(800, 400);
	image.create(WIDTH, HEIGHT);


	float maxdiv = 0;
	float avgdiv = 0;
	glm::vec2 avgdiv2 = { 0,0 };

	for (int y = 1; y < HEIGHT - 1; y++) {
		for (int x = 1; x < WIDTH - 1; x++) {
			float div = divergence(x, y, velocities);
			glm::vec2 div2 = divergence2(x, y, velocities);
			maxdiv = std::max(maxdiv, div2.x);
			avgdiv += div;
			avgdiv2 += div2;
		}
	}
	avgdiv2 /= (HEIGHT - 4) * (WIDTH - 4);
	std::cout << avgdiv2.x << " " << avgdiv2.y << " " << maxdiv << std::endl;

	for (int i = 0; i < 10; i++) {
		divergentVelocity();
	}


	maxdiv = 0;
	avgdiv = 0;
	avgdiv2 = { 0,0 };

	for (int y = 1; y < HEIGHT - 1; y++) {
		for (int x = 1; x < WIDTH - 1; x++) {
			float div = divergence(x, y, velocities);
			glm::vec2 div2 = divergence2(x, y, velocities);
			maxdiv = std::max(maxdiv, div2.x);
			avgdiv += div;
			avgdiv2 += div2;
		}
	}
	avgdiv2 /= (HEIGHT - 4) * (WIDTH - 4);
	std::cout << avgdiv2.x << " " << avgdiv2.y << " " << maxdiv << std::endl;

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		for (int y = HEIGHT / 2 - 50; y < HEIGHT / 2 + 50; y++) {
			for (int x = WIDTH / 2 - 50; x < WIDTH / 2 + 50; x++) {
				velocities[at(x,y)] = { 2,0 };
			}
		}

		convectback();



		for (int i = 0; i < dyes.size(); i++) {
			int x = i % WIDTH;
			int y = i / WIDTH;
			const sf::Color colour(
				/*(sf::Uint32)((velocities[i].r + 1) / 2 * 255),
				(sf::Uint32)((velocities[i].g + 1) / 2 * 255),
				0);*/
			(sf::Uint32)(dyes[i].r * 255),
			(sf::Uint32)(dyes[i].g * 255),
			(sf::Uint32)(dyes[i].b * 255));
			//std::cout << x << " " << y << " " << velocities[i].r << std::endl;
			image.setPixel(x, y, colour);
		}

		for (int i = 0; i < 10; i++) {
			divergentVelocity();
		}

		dyes = newDyes;
		for (auto& dye : newDyes) {
			dye = glm::vec3(0, 0, 0);
		}


		velocities = newVelocities;
		for (auto& vel : newVelocities) {
			vel = glm::vec2(0, 0);
		}

		fluidTexture.update(image);
		fluidSprite.setTexture(fluidTexture);

		window.clear();
		//window.draw(line, numVertex, sf::PrimitiveType::LinesStrip);
		window.draw(fluidSprite);
		window.display();
	}

	return 0;
}