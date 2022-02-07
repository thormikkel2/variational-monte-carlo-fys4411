#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <valarray>

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    
    double pos = 0;
    int dimensions = particles.at(0)->getNumberOfDimensions();
    
    for(int i = 0; i < particles.size(); i++){
        
        std::vector<double> currentPositions = particles.at(i)->getPosition();

        for(int j = 0; i < dimensions; j++){
            pos = pos + currentPositions.at(j)*currentPositions.at(j);
        }
    }    

    return exp(pos);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    //Need the sum of r_i^2 to calculate second derivative
    double sum_pos_squared = 0;
    int dimensions = particles.at(0)->getNumberOfDimensions();

    for(int i = 0; i < particles.size(); i++){
        
        std::vector<double> currentPositions = particles.at(i)->getPosition();

        for(int j = 0; i < dimensions; j++){
            sum_pos_squared = sum_pos_squared + currentPositions.at(j)*currentPositions.at(j);
        }
    }   

    //Return the analytical expression psi_T * 2a(2a*sum r_i^2 - 1)
    return exp(alpha*sum_pos_squared) * 2*alpha * (2*alpha*sum_pos_squared - 1));
}
