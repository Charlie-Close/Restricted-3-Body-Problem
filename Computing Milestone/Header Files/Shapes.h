//
//  Circle.h
//  OpenGlTest
//
//  Created by Charlie Close on 25/03/2023.
//

#ifndef Circle_h
#define Circle_h

#include<glad/glad.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>
#include<glm/gtx/rotate_vector.hpp>
#include<glm/gtx/vector_angle.hpp>
#include<glm/gtx/norm.hpp>
#include"Triangle.h"
#include"Mesh.h"


Mesh Circle(glm::vec3 position, float radius, int subDev);

#endif /* Circle_h */
