#include "assignment5.hpp"
#include "parametric_shapes.hpp"

#include "config.hpp"
#include "core/Bonobo.h"
#include "core/FPSCamera.h"
#include "core/helpers.hpp"
#include "core/node.hpp"
#include "core/ShaderProgramManager.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <imgui.h>
#include <tinyfiledialogs.h>

#include <clocale>
#include <stdexcept>

edaf80::Assignment5::Assignment5(WindowManager& windowManager) :
	mCamera(0.5f * glm::half_pi<float>(),
	        static_cast<float>(config::resolution_x) / static_cast<float>(config::resolution_y),
	        0.01f, 1000.0f),
	inputHandler(), mWindowManager(windowManager), window(nullptr)
{
	WindowManager::WindowDatum window_datum{ inputHandler, mCamera, config::resolution_x, config::resolution_y, 0, 0, 0, 0};

	window = mWindowManager.CreateGLFWWindow("EDAF80: Assignment 5", window_datum, config::msaa_rate);
	if (window == nullptr) {
		throw std::runtime_error("Failed to get a window: aborting!");
	}

	bonobo::init();
}

edaf80::Assignment5::~Assignment5()
{
	bonobo::deinit();
}

void
edaf80::Assignment5::run()
{
	// Set up the camera
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, 6.0f)); 
	mCamera.mMouseSensitivity = 0.003f;
	mCamera.mMovementSpeed = 3.0f; 
	auto camera_position = mCamera.mWorld.GetTranslation();
	glm::vec3 camera_start_angle = mCamera.mWorld.GetFront();

	////////////////////////////////////////////// Creating Shapes /////////////////////////////////////////

	// Cylinder shape
	float const cylinder_length = 1000.0f;
	float const cylinder_radius = 20.0f;
	auto cylinder_shape = parametric_shapes::createCylinder(cylinder_radius, cylinder_length, 500u, 500u);
	if (cylinder_shape.vao == 0u) {
		LogError("Failed to retrieve the mesh for the cylinder");
		return;
	}

	// Dragon shape
	auto dragon_shape = bonobo::loadObjects("/Users/feliciasegui/Desktop/Datorgrafik/CG_Labs/src/EDAF80/newdragon.obj")[0];

	// Torus shape
	float const torus_radius = 1.5f;
	auto torus_shape = parametric_shapes::createTorus(torus_radius, 0.15f, 30u, 30u);
	if (torus_shape.vao == 0u) {
		LogError("Failed to retrieve the mesh for the torus");
		return;
	}

	///////////////////////////////////////////////// Shaders ///////////////////////////////////////////////////
	
	ShaderProgramManager program_manager;

	// Dragon Phong Shader
	GLuint dragon_shader = 0u;
	program_manager.CreateAndRegisterProgram("DragonShader", 
											{ 	{ ShaderType::vertex, "EDAF80/phong_dragon.vert" } , 
												{ ShaderType::fragment, "EDAF80/phong_dragon.frag" } }, 
												dragon_shader);
	if (dragon_shader == 0u)
		LogError("Failed to load dragon shader");


	// Phong shader
	GLuint phong_shader = 0u;
	program_manager.CreateAndRegisterProgram("Phong", 
											{ 	{ ShaderType::vertex, "EDAF80/phong.vert" } , 
												{ ShaderType::fragment, "EDAF80/phong.frag" } }, 
												phong_shader);
	if (phong_shader == 0u)
		LogError("Failed to load phong shader");


	/////////////////////////////////////////////////////// Set uniforms ///////////////////////////////////////////////
	
	// Set dragon phong uniforms
	auto const light_position = glm::vec3(-2.0f, 4.0f, 2.0f);
	auto Ka = glm::vec3(0.0f);  
	auto Kd = glm::vec3(0.14f, 0.52f, 0.34f);
	auto Ks = glm::vec3(0.1f, 0.1f, 0.8f);
	auto illum = 10.0f;
	auto const set_dragon_uniforms = [&light_position, &camera_position, &Ka, &Kd, &Ks, &illum](GLuint program){
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "Ka"), 1, glm::value_ptr(Ka));
		glUniform3fv(glGetUniformLocation(program, "Kd"), 1, glm::value_ptr(Kd));
		glUniform1f(glGetUniformLocation(program, "illum"), illum);
	};

	// Set Phong Uniforms for Tori
	bool use_normal_mapping = true;
	auto ambient = glm::vec3(0.1f, 0.1f, 0.1f);
	auto diffuse = glm::vec3(1.0f, 1.0f, 1.0f);
	auto specular = glm::vec3(1.0f, 1.0f, 1.0f);
	auto shininess = 10.0f;
	auto const phong_set_uniforms_tori = [&use_normal_mapping,&light_position,&camera_position,&ambient,&diffuse,&specular,&shininess](GLuint program){
		glUniform1i(glGetUniformLocation(program, "use_normal_mapping"), use_normal_mapping ? 1 : 0);
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "ambient"), 1, glm::value_ptr(ambient));
		glUniform3fv(glGetUniformLocation(program, "diffuse"), 1, glm::value_ptr(diffuse));
		glUniform3fv(glGetUniformLocation(program, "specular"), 1, glm::value_ptr(specular));
		glUniform1f(glGetUniformLocation(program, "shininess"), shininess);
	};

	// Set Phong Uniforms for front torus
	ambient = glm::vec3(0.1f, 0.1f, 0.1f);
	auto const phong_set_uniforms_front_torus = [&use_normal_mapping,&light_position,&camera_position,&ambient,&diffuse,&specular,&shininess](GLuint program){
		glUniform1i(glGetUniformLocation(program, "use_normal_mapping"), use_normal_mapping ? 1 : 0);
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "ambient"), 1, glm::value_ptr(ambient));
		glUniform3fv(glGetUniformLocation(program, "diffuse"), 1, glm::value_ptr(diffuse));
		glUniform3fv(glGetUniformLocation(program, "specular"), 1, glm::value_ptr(specular));
		glUniform1f(glGetUniformLocation(program, "shininess"), shininess);
	};

	// Set cylinder Phong uniforms 
	auto ambientcyl = glm::vec3(0.7f, 0.1f, 0.1f);
	bool use_normal_mapping_cyl = false;
	auto const phong_set_uniforms_cylinder = [&use_normal_mapping_cyl,&light_position,&camera_position,&ambientcyl,&diffuse,&specular,&shininess](GLuint program){
		glUniform1i(glGetUniformLocation(program, "use_normal_mapping"), use_normal_mapping_cyl ? 1 : 0);
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "ambient"), 1, glm::value_ptr(ambientcyl));
		glUniform3fv(glGetUniformLocation(program, "diffuse"), 1, glm::value_ptr(diffuse));
		glUniform3fv(glGetUniformLocation(program, "specular"), 1, glm::value_ptr(specular));
		glUniform1f(glGetUniformLocation(program, "shininess"), shininess);
	};

	///////////////////////////////////////////////// Textures //////////////////////////////////////////////////
	
	// Phong Lava textures 
	auto my_diffuse_tex_id = bonobo::loadTexture2D("/Users/feliciasegui/Desktop/Datorgrafik/CG_Labs/res/textures/lava_diff.png");
	auto my_specular_map_id = bonobo::loadTexture2D("/Users/feliciasegui/Desktop/Datorgrafik/CG_Labs/res/textures/lava_spec.png");
	auto my_normal_map_id = bonobo::loadTexture2D("/Users/feliciasegui/Desktop/Datorgrafik/CG_Labs/res/textures/lava_normal.png");

	///////////////////////////////////////////// Creating Nodes ///////////////////////////////////////////////////

	// Cylinder Locations
	int const nbr_cylinders = 3;
	std::array<glm::vec3, nbr_cylinders> cylinder_locations = {
		glm::vec3(0.0f,  0.0f,  0.0f),
		glm::vec3(0.0f,  0.0f,  - 1.0f * cylinder_length),
		glm::vec3(0.0f,  0.0f,  - 2.0f * cylinder_length)
	};
	std::array<glm::vec3, nbr_cylinders> cylinder_start_locations = cylinder_locations;

	// Node for first cylinder
	std::array<Node, nbr_cylinders> cylinders;
	
	auto& cylinder = cylinders[0];
	cylinder.set_geometry(cylinder_shape);
	cylinder.set_program(&phong_shader, phong_set_uniforms_cylinder);
	cylinder.get_transform().SetTranslate(cylinder_locations[0]);
	cylinder.get_transform().SetRotate(glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	cylinder.add_texture("diffuse_text", my_diffuse_tex_id, GL_TEXTURE_2D);
	cylinder.add_texture("specular_map", my_specular_map_id, GL_TEXTURE_2D);
	cylinder.add_texture("normal_map", my_normal_map_id, GL_TEXTURE_2D);

	glm::vec3 previous_cylinder_position = cylinder.get_transform().GetTranslation();

	// Nodes for the rest of the cylinders
	for (std::size_t i = 1; i < nbr_cylinders; ++i) {
		auto& cylinder = cylinders[i];

		cylinder.set_geometry(cylinder_shape);
		cylinder.set_program(&phong_shader, phong_set_uniforms_cylinder);

		cylinder.add_texture("diffuse_text", my_diffuse_tex_id, GL_TEXTURE_2D);
		cylinder.add_texture("specular_map", my_specular_map_id, GL_TEXTURE_2D);
		cylinder.add_texture("normal_map", my_normal_map_id, GL_TEXTURE_2D);

		glm::vec3 new_position = previous_cylinder_position + glm::vec3(0.0f, 0.0f, - cylinder_length);

		cylinder_locations[i] = new_position;
		cylinder.get_transform().SetTranslate(new_position);
		cylinder.get_transform().SetRotate(glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));

		previous_cylinder_position = cylinder.get_transform().GetTranslation();
	}

	// Node dragon
	Node dragon;
	dragon.set_geometry(dragon_shape);
	dragon.set_program(&dragon_shader, set_dragon_uniforms);
	glm::vec3 dragon_position = glm::vec3(0.0f, 0.0f, 0.0f);
	dragon.get_transform().SetTranslate(dragon_position);
	

	// Node tori
	int const nbr_tori = 5;
	std::array<glm::vec3, nbr_tori> tori_locations = {
		glm::vec3(0.0f,  0.0f,  0.0f),
		glm::vec3(1.0f,  0.0f,  1.0f),
		glm::vec3(2.0f,  0.0f,  2.0f),
		glm::vec3(3.0f,  0.0f,  3.0f),
		glm::vec3(4.0f,  1.0f,  3.0f),
	};
	std::array<glm::vec3, nbr_tori> tori_start_locations = tori_locations;

	std::array<Node, nbr_tori> tori;

	auto& torus = tori[0];
	torus.set_geometry(torus_shape);
	torus.set_program(&phong_shader, phong_set_uniforms_tori);
	torus.get_transform().SetTranslate(glm::vec3(0.0f, 0.0f, 0.0f));
	torus.get_transform().SetRotateX(glm::half_pi<float>());
	

	srand(time(0));
	float radious = 12;
	glm::vec3 previous_torus_position = torus.get_transform().GetTranslation();

	for (std::size_t i = 1; i < nbr_tori; ++i) {
		auto& torus = tori[i];
		torus.set_geometry(torus_shape);
		torus.set_program(&phong_shader, phong_set_uniforms_front_torus);
		torus.add_texture("diffuse_text", my_diffuse_tex_id, GL_TEXTURE_2D);
		torus.add_texture("specular_map", my_specular_map_id, GL_TEXTURE_2D);
		torus.add_texture("normal_map", my_normal_map_id, GL_TEXTURE_2D);

		float theta = glm::two_pi<float>() * static_cast<float>(rand() % 101) / 100;
		float phi = glm::pi<float>() * static_cast<float>(rand() % 51 +50) / 600;
		glm::vec3 new_position = previous_torus_position + glm::vec3(radious * sin(phi) * cos(theta), radious * sin(phi) * sin(theta), -radious * cos(phi));
		
		tori_locations[i] = new_position;
		torus.get_transform().SetTranslate(new_position);
		torus.get_transform().SetRotateX(glm::half_pi<float>());

		previous_torus_position = torus.get_transform().GetTranslation();
	}

	/////////////////////////////////////// Setup values ///////////////////////////////////////////
	bool gameover = false;
	float ellapsed_time_s = 0.0f;

	// Score 
	int score = 0;
	int highscore = 0;
	bool new_possibility_to_score = true;
	
	glClearDepthf(1.0f);
	glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	glEnable(GL_DEPTH_TEST);

	auto lastTime = std::chrono::high_resolution_clock::now();

	auto cull_mode = bonobo::cull_mode_t::disabled;
	auto polygon_mode = bonobo::polygon_mode_t::fill;
	bool show_logs = false;
	bool show_gui = true;
	bool shader_reload_failed = false;
	bool show_basis = false;
	float basis_thickness_scale = 1.0f;
	float basis_length_scale = 1.0f;

	changeCullMode(cull_mode);

	int current_state = PLAY_GAME; 

	while (!glfwWindowShouldClose(window)) {
		switch (current_state)
		{
		//////////////////// NEW GAME //////////////////////
		case 1: {

			// Camera
			mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, 6.0f)); 
			mCamera.mWorld.LookTowards(camera_start_angle);
			auto camera_position = mCamera.mWorld.GetTranslation();

			// Cylinders
			for (std::size_t i = 0; i < nbr_cylinders; ++i) {
				auto& cylinder = cylinders[i];
				glm::vec3 new_position = cylinder_start_locations[i];
				cylinder_locations[i] = new_position;
				cylinder.get_transform().SetTranslate(new_position);
				cylinder.get_transform().SetRotate(glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
			}

			// Dragon
			dragon.get_transform().SetTranslate(glm::vec3(0.0f, 0.0f, 0.0f));
			inputHandler.SetVelocity(glm::vec3(0.0f));

			// First Torus
			auto& torus = tori[0];
			glm::vec3 new_position = tori_start_locations[0];
			torus.get_transform().SetTranslate(new_position);
			torus.get_transform().SetRotateX(glm::half_pi<float>());
			previous_torus_position = new_position;
			tori_locations[0] = new_position;
			
			// The rest of the tori
			for (std::size_t i = 1; i < nbr_tori; ++i) {
				auto& torus = tori[i];
				float theta = glm::two_pi<float>() * static_cast<float>(rand() % 101) / 100;
				float phi = glm::pi<float>() * static_cast<float>(rand() % 51 +50) / 600;
				glm::vec3 new_position = previous_torus_position + glm::vec3(radious * sin(phi) * cos(theta), radious * sin(phi) * sin(theta), -radious * cos(phi));
				tori_locations[i] = new_position;
				torus.get_transform().SetTranslate(new_position);
				torus.get_transform().SetRotateX(glm::half_pi<float>());
				previous_torus_position = new_position;
			}
			
			// Score 
			score = 0;
			new_possibility_to_score = true;
			ellapsed_time_s = 0.0f;
			gameover = false;
			current_state = PLAY_GAME;
			break;

		} 
		////////////////////// PLAY GAME //////////////////////////////
		case 2: {
			auto const nowTime = std::chrono::high_resolution_clock::now();
			auto const deltaTimeUs = std::chrono::duration_cast<std::chrono::microseconds>(nowTime - lastTime);
			lastTime = nowTime;
			ellapsed_time_s += std::chrono::duration<float>(deltaTimeUs).count();

			auto& io = ImGui::GetIO();
			inputHandler.SetUICapture(io.WantCaptureMouse, io.WantCaptureKeyboard);

			glfwPollEvents();
			inputHandler.Advance();
			mCamera.Update(deltaTimeUs, inputHandler);
			camera_position = mCamera.mWorld.GetTranslation();

			////////////////////// Create a new cylinder when the first is behind us ///////////////////////
			if(camera_position.z < (cylinder_locations[0].z-cylinder_length)){
				auto cylinder = cylinders[0];
			
				glm::vec3 new_position = cylinder_locations[nbr_cylinders-1] + glm::vec3(0.0f, 0.0f, -cylinder_length);
				cylinder.get_transform().SetTranslate(new_position);
				cylinder.get_transform().SetRotate(glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));

				for (int i = 0; i < nbr_cylinders - 1; i++){
					cylinder_locations[i] = cylinder_locations[i+1];
					cylinders[i] = cylinders[i+1];
				} 
				cylinder_locations[nbr_cylinders-1] = new_position;
				cylinders[nbr_cylinders-1] = cylinder;
			}

			///////////////////// Create a new torus when the first is behind us //////////////////////
			if (glm::dot(mCamera.mWorld.GetFront(), tori_locations[0] - camera_position) < 0) {
			
				new_possibility_to_score = true;

				previous_torus_position = tori_locations[nbr_tori - 1];
				float theta = glm::two_pi<float>() * static_cast<float>(rand() % 101) / 100;
				float phi = glm::pi<float>() * static_cast<float>(rand() % 51 + 50) / 600;
				glm::vec3 new_position = previous_torus_position + glm::vec3(radious * sin(phi) * cos(theta), radious * sin(phi) * sin(theta), -radious * cos(phi));

				while(glm::length(glm::vec2(new_position.x, new_position.y)) > (cylinder_radius - torus_radius)){
					theta = glm::two_pi<float>() * static_cast<float>(rand() % 101) / 100;
					phi = glm::pi<float>() * static_cast<float>(rand() % 51 + 50) / 600;
					new_position = previous_torus_position + glm::vec3(radious * sin(phi) * cos(theta), radious * sin(phi) * sin(theta), -radious * cos(phi));
				}

				for (int i = 0; i < nbr_tori-1; ++i) {
					auto& torus = tori[i]; // 
					tori_locations[i] = tori_locations[i+1];
					torus.get_transform().SetTranslate(tori_locations[i]);
					torus.get_transform().SetRotateX(glm::half_pi<float>());

				}
				auto& torus = tori[nbr_tori-1];
				tori_locations[nbr_tori-1] = new_position;
				torus.get_transform().SetTranslate(tori_locations[nbr_tori-1]);
				torus.get_transform().SetRotateX(glm::half_pi<float>());
			}

			////////////////////////////// Translate dragon /////////////////////////////////////

			// Translate the player object to always being in front of the camera
			if (!gameover) dragon_position =  camera_position + mCamera.mWorld.GetFront() * 3.0f - glm::vec3(0.0f, 0.11f, 0.0f);
	
			dragon.get_transform().SetTranslate(dragon_position);
		
			// Rotate the dragon so that it always looks away from the camera
			dragon.get_transform().LookAt(camera_position);
			dragon.get_transform().Rotate(glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));

			// Create parent matrix from the scale and translation matrix
			glm::vec3 const dragon_scale{ 0.0008f };
			glm::mat4 scale = glm::scale(glm::mat4(1.0f), dragon_scale);
			glm::mat4 trans_matrix = dragon.get_transform().GetTranslationMatrix();
			glm::mat4 world = trans_matrix*scale ;
		
			////////////////////////// Check if we are inside a torus -> score /////////////////////////

			glm::vec3 torus_loc = tori_locations[0];
			float pcp_pd = glm::length(torus_loc-dragon_position);
			if (pcp_pd <= torus_radius && new_possibility_to_score){
				score +=1;
				if (highscore < score ) highscore = score;
				new_possibility_to_score = false;
			}

			////////////////////// Check if we are outside our world -> game over ///////////////////////////////
		
			if (glm::length(glm::vec2(dragon_position.x, dragon_position.y)) > cylinder_radius) gameover = true;

			/////////////////////////// Input handler ///////////////////////////////////////
			if (inputHandler.GetKeycodeState(GLFW_KEY_R) & JUST_PRESSED) {
				shader_reload_failed = !program_manager.ReloadAllPrograms();
				if (shader_reload_failed)
					tinyfd_notifyPopup("Shader Program Reload Error",
				                   "An error occurred while reloading shader programs; see the logs for details.\n"
				                   "Rendering is suspended until the issue is solved. Once fixed, just reload the shaders again.",
				                   "error");
			}
			if (inputHandler.GetKeycodeState(GLFW_KEY_F3) & JUST_RELEASED)
				show_logs = !show_logs;
			if (inputHandler.GetKeycodeState(GLFW_KEY_F2) & JUST_RELEASED)
				show_gui = !show_gui;
			if (inputHandler.GetKeycodeState(GLFW_KEY_F11) & JUST_RELEASED)
				mWindowManager.ToggleFullscreenStatusForWindow(window);

			int framebuffer_width, framebuffer_height;
			glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);
			glViewport(0, 0, framebuffer_width, framebuffer_height);

			mWindowManager.NewImGuiFrame();

			glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
			bonobo::changePolygonMode(polygon_mode);

			////////////////////////////////// Render //////////////////////////////////////////
			if (!shader_reload_failed) {

				cylinder.render(mCamera.GetWorldToClipMatrix());
				dragon.render(mCamera.GetWorldToClipMatrix(), world);
				for (std::size_t i = 0; i < nbr_tori; ++i) {
					auto& torus = tori[i];
					torus.render(mCamera.GetWorldToClipMatrix());
				}
				for (std::size_t i = 0; i < nbr_cylinders; i++){
					auto& cylinder = cylinders[i];
					cylinder.render(mCamera.GetWorldToClipMatrix());
				}

			}

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

			/////////////////////////////////// ImGui /////////////////////////////////////////

			bool opened = ImGui::Begin("The Eternal Tornment of the Dragon", nullptr, ImGuiWindowFlags_None);
			if (opened) {
				ImGui::Text("Time: %f", ellapsed_time_s);
				ImGui::Text("Session Highscore: %d", highscore);
				ImGui::TextColored(ImVec4(1,1,0,1), "Score: %d", score);
				//ImGui::Text("X pos: %f\nY pos: %f\nZ pos: %f", dragon_position.x, dragon_position.y, dragon_position.z );
				if (gameover){
					ImGui::Text("Game over, want to play again? \nY = yes, N = no\n");
					if (inputHandler.GetKeycodeState(GLFW_KEY_Y) & JUST_RELEASED) {
						gameover = false;
						current_state = NEW_GAME;
					} else if (inputHandler.GetKeycodeState(GLFW_KEY_N) & JUST_RELEASED) {
						current_state = END_GAME;
					}
				}
			
			}
			ImGui::End();

			mWindowManager.RenderImGuiFrame(show_gui);

			glfwSwapBuffers(window);

			break;


		} 
		////////////////////////////////////// END GAME //////////////////////////////////////////////
		case 3: {
			printf("-----------------------------\n        Game Over\n    Session High Score: %d\n-----------------------------\n", highscore);
			glfwSetWindowShouldClose(window, GLFW_TRUE);
			break;
		}

		}


	}


}

int main()
{
	std::setlocale(LC_ALL, "");

	Bonobo framework;

	try {
		edaf80::Assignment5 assignment5(framework.GetWindowManager());
		assignment5.run();	
	} catch (std::runtime_error const& e) {
		LogError(e.what());
	}
}
