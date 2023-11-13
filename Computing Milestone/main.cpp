//
//  main.cpp
//  Computing Milestone
//
//  Created by Charlie Close on 28/10/2023.
//

#include <iostream>
#include <thread>
#include "Mesh.h"
#include "Shapes.h"
#include <fstream>


const long double M_e = 5.9742 * pow(10, 24);
const long double M_m = 7.35 * pow(10, 22);
const long double G = 6.6726 * pow(10, -11);
const long double D = 3.844 * pow(10, 5) * pow(10, 3); //meters
const long double r_e = D * M_m / (M_m + M_e);
const long double r_m = D * M_e / (M_m + M_e);
const long double GM_e = G * M_e;
const long double GM_m = G * M_m;
const unsigned int width = 1200;
const unsigned int height = 800;

const long double a_m = GM_e / pow(r_m, 2);
const long double dt = 0.0005; //seconds
const long double dt2_2 = pow(dt, 2) / 2;
const long double dt_2 = dt / 2;
const long double dt_6 = dt / 6;
const long double time_to_simulate = 20 * 28 * 24 * 60 * 60;
const int total_steps = round(time_to_simulate / dt);
const long double dtheta = pow(a_m / r_m, 0.5) * dt;
const long double dtheta_2 = dtheta / 2;

const long double initL2 = D * pow(M_m / (3 * M_e), double(1)/double(3));
const int cpus = 10;

const float shrink_factor = 100000;


struct positions {
    long double e[3];
    long double m[3];
};

struct accelerationType {
    long double ac[3];
};

positions findPositions(long double theta) {
    long double cos_t = cos(theta);
    long double sin_t = sin(theta);

    positions poses;
    poses.e[0] = - cos_t * r_e;
    poses.e[1] = - sin_t * r_e;
    poses.e[2] = 0;
    
    poses.m[0] = cos_t * r_m;
    poses.m[1] = sin_t * r_m;
    poses.m[2] = 0;

    return poses;
}

double distToL2(long double P_m[3], long double P[3], long double L2) {
    double dist_m = pow(pow(P_m[0] - P[0], 2) + pow(P_m[1] - P[1], 2) + pow(P_m[2] - P[2], 2), 0.5);
    return dist_m - L2;
}

accelerationType calculateAcceleration(long double position[3], positions poses) {
    long double v_e[3] = { position[0] - poses.e[0], position[1] - poses.e[1], position[2] - poses.e[2] };
    long double v_m[3] = { position[0] - poses.m[0], position[1] - poses.m[1], position[2] - poses.m[2] };
    
    long double d_e2 = pow(v_e[0], 2) + pow(v_e[1], 2) + pow(v_e[2], 2);
    long double d_m2 = pow(v_m[0], 2) + pow(v_m[1], 2) + pow(v_m[2], 2);
    
    long double e_mult = GM_e / pow(d_e2, 1.5);
    long double m_mult = GM_m / pow(d_m2, 1.5);
    
    long double acceleration[3] = { - e_mult * v_e[0] - m_mult * v_m[0], - e_mult * v_e[1] - m_mult * v_m[1], - e_mult * v_e[2] - m_mult * v_m[2] };
    accelerationType toRet;
    toRet.ac[0] = acceleration[0];
    toRet.ac[1] = acceleration[1];
    toRet.ac[2] = acceleration[2];
    return toRet;
}

int runTaylor(long double L2, int steps, long double *val) {
    long double position[3] = {r_m + L2, 0, 0};
    long double velocity[3] = {0, position[0] * dtheta / dt, 0};
    
    int step = 0;
    long double theta = 0;
    
    while (step < steps) {
        step += 1;
        positions poses = findPositions(theta);
        
        accelerationType acc = calculateAcceleration(position, poses);
        
        position[0] += dt * velocity[0] + dt2_2 * acc.ac[0];
        position[1] += dt * velocity[1] + dt2_2 * acc.ac[1];
        position[2] += dt * velocity[2] + dt2_2 * acc.ac[2];
        
        velocity[0] += dt * acc.ac[0];
        velocity[1] += dt * acc.ac[1];
        velocity[2] += dt * acc.ac[2];
        
        theta += dtheta;
    }
    *val = distToL2(findPositions(theta).m, position, L2);
    return 0;
}

int runTaylor_Plot(long double L2, int steps) {
    long double position[3] = {r_m + L2, 0, 0};
    long double velocity[3] = {0, position[0] * dtheta / dt, 0};
    
    int step = 0;
    long double theta = 0;
    
    std::vector<float> dataPoints;
    
    const int skips = steps / 5000;
    
    while (step < steps) {
        step += 1;
        positions poses = findPositions(theta);
        
        if (step % skips == 0) {
            dataPoints.push_back(float(distToL2(poses.m, position, L2)));
        }
        
        accelerationType acc = calculateAcceleration(position, poses);
        
        position[0] += dt * velocity[0] + dt2_2 * acc.ac[0];
        position[1] += dt * velocity[1] + dt2_2 * acc.ac[1];
        position[2] += dt * velocity[2] + dt2_2 * acc.ac[2];
        
        velocity[0] += dt * acc.ac[0];
        velocity[1] += dt * acc.ac[1];
        velocity[2] += dt * acc.ac[2];
        
        theta += dtheta;
    }
    
    std::cout << "Writing File!" << std::endl;
    std::string data = "";
    
    for (int i = 0; i < dataPoints.size(); i++) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(30) << dataPoints[i];
        std::string s = stream.str();
        
        data += s + ", " + std::to_string(dt * (double(i) * double(skips + 1))) + "\n";
    }
    
    std::ofstream MyFile("/Users/charlieclose/Documents/Programs/C++/Computing Milestone/Computing Milestone/TS.txt");
    MyFile << data;
    MyFile.close();
    return 0;
}

int runRK4(long double L2, int steps, long double *val) {
    long double position[3] = {r_m + L2, 0, 0};
    long double velocity[3] = {0, position[0] * dtheta / dt, 0};
    
    int step = 0;
    long double theta = 0;
    
    while (step < steps) {
        step += 1;
    
        const positions pos1 = findPositions(theta + dtheta_2);
        long double z1[3] = { position[0] + dt_2 * velocity[0], position[1] + dt_2 * velocity[1], position[2] + dt_2 * velocity[2] };
        accelerationType acc = calculateAcceleration(position, pos1);
        long double dz1[3] = { velocity[0] + dt_2 * acc.ac[0], velocity[1] + dt_2 * acc.ac[1], velocity[2] + dt_2 * acc.ac[2] };
        
        long double z2[3] = { position[0] + dt_2 * dz1[0], position[1] + dt_2 * dz1[1], position[2] + dt_2 * dz1[2] };
        accelerationType acc1 = calculateAcceleration(z1, pos1);
        long double dz2[3] = { velocity[0] + dt_2 * acc1.ac[0], velocity[1] + dt_2 * acc1.ac[1], velocity[2] + dt_2 * acc1.ac[2] };
        
        const positions pos2 = findPositions(theta + dtheta);
        long double z3[3] = { position[0] + dt * dz2[0], position[1] + dt * dz2[1], position[2] + dt * dz2[2] };
        accelerationType acc2 = calculateAcceleration(z2, pos2);
        long double dz3[3] = { velocity[0] + dt * acc2.ac[0], velocity[1] + dt * acc2.ac[1], velocity[2] + dt * acc2.ac[2] };
        
        accelerationType acc3 = calculateAcceleration(z3, pos2);
        
        position[0] += dt_6 * (velocity[0] + 2 * dz1[0] + 2 * dz2[0] + dz3[0]);
        position[1] += dt_6 * (velocity[1] + 2 * dz1[1] + 2 * dz2[1] + dz3[1]);
        position[2] += dt_6 * (velocity[2] + 2 * dz1[2] + 2 * dz2[2] + dz3[2]);
        
        velocity[0] += dt_6 * (acc.ac[0] + 2 * acc1.ac[0] + 2 * acc2.ac[0] + acc3.ac[0]);
        velocity[1] += dt_6 * (acc.ac[1] + 2 * acc1.ac[1] + 2 * acc2.ac[1] + acc3.ac[1]);
        velocity[2] += dt_6 * (acc.ac[2] + 2 * acc1.ac[2] + 2 * acc2.ac[2] + acc3.ac[2]);
        
        theta += dtheta;
    }
    *val = distToL2(findPositions(theta).m, position, L2);
    return 0;
}

int runRK4_Plot(long double L2, int steps) {
    long double position[3] = {r_m + L2, 0, 0};
    long double velocity[3] = {0, position[0] * dtheta / dt, 0};
    
    int step = 0;
    long double theta = 0;
    
    std::vector<float> dataPoints;
    
    const int skips = steps / 5000;
    
    while (step < steps) {
        step += 1;
    
        positions pos1 = findPositions(theta + dtheta_2);
        long double z1[3] = { position[0] + dt_2 * velocity[0], position[1] + dt_2 * velocity[1], position[2] + dt_2 * velocity[2] };
        accelerationType acc = calculateAcceleration(position, pos1);
        long double dz1[3] = { velocity[0] + dt_2 * acc.ac[0], velocity[1] + dt_2 * acc.ac[1], velocity[2] + dt_2 * acc.ac[2] };
        
        long double z2[3] = { position[0] + dt_2 * dz1[0], position[1] + dt_2 * dz1[1], position[2] + dt_2 * dz1[2] };
        accelerationType acc1 = calculateAcceleration(z1, pos1);
        long double dz2[3] = { velocity[0] + dt_2 * acc1.ac[0], velocity[1] + dt_2 * acc1.ac[1], velocity[2] + dt_2 * acc1.ac[2] };
        
        positions pos2 = findPositions(theta + dtheta);
        long double z3[3] = { position[0] + dt * dz2[0], position[1] + dt * dz2[1], position[2] + dt * dz2[2] };
        accelerationType acc2 = calculateAcceleration(z2, pos2);
        long double dz3[3] = { velocity[0] + dt * acc2.ac[0], velocity[1] + dt * acc2.ac[1], velocity[2] + dt * acc2.ac[2] };
        
        accelerationType acc3 = calculateAcceleration(z3, pos2);
        
        position[0] += dt_6 * (velocity[0] + 2 * dz1[0] + 2 * dz2[0] + dz3[0]);
        position[1] += dt_6 * (velocity[1] + 2 * dz1[1] + 2 * dz2[1] + dz3[1]);
        position[2] += dt_6 * (velocity[2] + 2 * dz1[2] + 2 * dz2[2] + dz3[2]);
        
        velocity[0] += dt_6 * (acc.ac[0] + 2 * acc1.ac[0] + 2 * acc2.ac[0] + acc3.ac[0]);
        velocity[1] += dt_6 * (acc.ac[1] + 2 * acc1.ac[1] + 2 * acc2.ac[1] + acc3.ac[1]);
        velocity[2] += dt_6 * (acc.ac[2] + 2 * acc1.ac[2] + 2 * acc2.ac[2] + acc3.ac[2]);
        
        theta += dtheta;
        
        if (step % skips == 0) {
            dataPoints.push_back(float(distToL2(pos2.m, position, L2)));
        }
    }
    
    std::cout << "Writing File!" << std::endl;
    std::string data = "";
    
    for (int i = 0; i < dataPoints.size(); i++) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(30) << dataPoints[i];
        std::string s = stream.str();
        
        data += s + ", " + std::to_string(dt * (i * skips + 1)) + "\n";
    }
    
    std::ofstream MyFile("/Users/charlieclose/Documents/Programs/C++/Computing Milestone/Computing Milestone/RK4.txt");
    MyFile << data;
    MyFile.close();
    
    return 0;
}

int main() {
    int steps = 200000000;
    std::cout << std::setprecision(150) << (50.f * float(steps) * dt) / (60.f * 60.f * 24.f) << std::endl;
    
    long double lowerBound = initL2 * 0.9;
    long double upperBound = initL2 * 1.1;
    
    for (int i = 0; i < 20; i++) {
        long double boundStep = (upperBound - lowerBound) / (cpus + 1);
        long double bounds[cpus];
        long double values[cpus];
        std::thread threads[cpus];
        for (int cpu = 0; cpu < cpus; cpu ++) {
            bounds[cpu] = lowerBound + (cpu + 1) * boundStep;
            threads[cpu] = std::thread(runTaylor, bounds[cpu], steps, &values[cpu]);
        }
        bool found = false;
        for (int cpu = 0; cpu < cpus; cpu ++) {
            threads[cpu].join();
            if (!found) {
                if (values[cpu] < 0) {
                    lowerBound = bounds[cpu];
                } else {
                    found = true;
                    upperBound = bounds[cpu];
                }
            }
        }
        std::cout << i << std::endl;
    }
    
    const long double average = (upperBound + lowerBound) / 2;
    
    runTaylor_Plot(average, steps * 50);
    
    std::ofstream MyFile("/Users/charlieclose/Documents/Programs/C++/Computing Milestone/Computing Milestone/TSValues.txt");
    std::stringstream stream;
    stream << std::fixed << std::setprecision(60) << (upperBound + lowerBound) / 2;
    std::string value = stream.str();
    
    std::stringstream stream2;
    stream2 << std::fixed << std::setprecision(60) << (upperBound - lowerBound) / 2;
    std::string error = stream2.str();
    
    std::string data = "L2 Point: " + value + "\nError: " + error;
    MyFile << data;
//    runWithGraphics((lowerBound + upperBound) / 2);
    return 0;
}
