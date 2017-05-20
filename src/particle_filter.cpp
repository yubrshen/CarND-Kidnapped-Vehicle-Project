/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <map>

#include "particle_filter.h"

using namespace std; // required for access default_random_engine
static default_random_engine gen(4211692);

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// DONE: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  num_particles = 20; // give it a try to see if it's fast enough, and accurate enough

	double std_x, std_y, std_theta; // Standard deviations for x, y, and psi

  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];

  // Gaussian random number generators with given mean and std
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

  for (uint i = 0; i < num_particles; ++i) {
    Particle particle; // statically allocated in the scope of the ParticleFilter life span, should be sufficient of life-span.
    particle.id = i;

    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;

    particles.push_back(particle);
    weights.push_back(particle.weight);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// DONE: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	double std_x, std_y, std_theta; // Standard deviations for x, y, and psi
  double vel_yaw = velocity/yaw_rate;
  double yaw_dt = yaw_rate*delta_t;
  double vel_dt = velocity*delta_t;

  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2];
  // Gaussian random number generators with given mean and std
  normal_distribution<double> dist_x(0.0, std_x); // the mean should be 0.0 for noise
  normal_distribution<double> dist_y(0.0, std_y);
  normal_distribution<double> dist_theta(0.0, std_theta);

  for (uint i = 0; i < num_particles; ++i) {
    // motion model:
    if (0.0001 < fabs(yaw_rate)) {
      const double new_theta = particles[i].theta + yaw_dt;
      particles[i].x += vel_yaw * ( sin(new_theta) - sin(particles[i].theta));
      particles[i].y += vel_yaw * (-cos(new_theta) + cos(particles[i].theta));
      particles[i].theta = new_theta;
    } else {
      particles[i].x += vel_dt * cos(particles[i].theta);
      particles[i].y += vel_dt * sin(particles[i].theta);
      // no change to particles[i].theta
    }
    // add random noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

// The following is no longer needed. Its functionality is absorbed in updateWeights.
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// DONE: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  // Implementation outline:
  // For each particle:
  // For each observation:
  // transform the observation to global coordinate system
  // find the closest matched landmark for the observation
  // only consider those landmarks within the sensor range
  // compute the probability of the observation for the landmark
  // compute the overall weight-importance for all the observations for the particle, normalize the weights

  Map::single_landmark_s observation_global_coord;
  uint closest_landmark_index = -1;
  double dist_tmp = 0;
  double smallest_distance = 0;
  Map::single_landmark_s matched_landmark;
  const double sigma_x_squared = std_landmark[0] * std_landmark[0];
  const double sigma_y_squared = std_landmark[1] * std_landmark[1];
  //const double k = 2 * M_PI * std_landmark[0] * std_landmark[1];

  double dx = 0.0;
  double dy = 0.0;
  double weights_total = 0.0;

  for (uint i = 0; i < num_particles; ++i) {
    double weight_superscript = 0.0;
    const double sin_theta = sin(particles[i].theta);
    const double cos_theta = cos(particles[i].theta);

    for (uint j = 0; j < observations.size(); ++j) {
      // transform observation from local coordinate to global one
      observation_global_coord.x_f = particles[i].x + (observations[j].x * cos_theta) - (observations[j].y * sin_theta);
      observation_global_coord.y_f = particles[i].y + (observations[j].x * sin_theta) + (observations[j].y * cos_theta);

      closest_landmark_index = -1;
      dist_tmp = 0;
      smallest_distance = 0;
      for (uint l = 0; l < map_landmarks.landmark_list.size(); ++l) {
        dist_tmp = dist(observation_global_coord.x_f, observation_global_coord.y_f,
                        map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f);
        if (closest_landmark_index == -1 || dist_tmp < smallest_distance) {
          smallest_distance = dist_tmp;
          closest_landmark_index = l;
        }
      }

      if (smallest_distance < sensor_range) {
        matched_landmark = map_landmarks.landmark_list[closest_landmark_index];
        // compute the probability of the observation in relation to the matched landmark, by Gaussian distribution,
        // and accumulate it to the multiplication of probabilities
        dx = observation_global_coord.x_f - matched_landmark.x_f;
        dy = observation_global_coord.y_f - matched_landmark.y_f;
        weight_superscript += dx*dx/sigma_x_squared + dy*dy/sigma_y_squared;
        // Optimization: the multiplication of exponentials can be computed by first add up the superscripts.
        // then apply the summed subscript with exponent at the end below.
      }
    }
    particles[i].weight = exp(-0.5*weight_superscript);
    // Optimization: no need to divide the above by 2 * M_PI * std_landmark[0] * std_landmark[1],
    // a constant for Gaussian distribution
    // as the weights will be normalized by division of the total sum of the weights.
    weights_total += particles[i].weight;
  }
  // normalize weights to sum 1
  for (uint i = 0; i < num_particles; ++i) {
    particles[i].weight = particles[i].weight / weights_total;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
	// DONE: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::vector<Particle> new_particles;

  std::discrete_distribution<int> particle_distribution(weights.begin(), weights.end());
  std::map<int, int> m;
  uint selected_id = 0;
  for (uint i = 0; i < num_particles; ++i) {
    selected_id = particle_distribution(gen);
    ++m[selected_id];
    new_particles.push_back(particles[selected_id]);
  }
  cout << "picked sample number: " << m.size() << endl;
  for (auto p : m) {
    cout << p.first << " selected with probability " << float(p.second)/num_particles << " the particle probability: " << particles[p.first].weight << endl;
  };
  particles = new_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
