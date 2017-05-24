#ifndef _RNG_H
#define _RNG_H

#include <random>

static const float min_size = 0.001;

class rng_distribution {
  public:
	virtual ~rng_distribution() {}
	virtual float operator()() = 0;

	void seed(int s) { rng.seed(s); }

  protected:
	std::default_random_engine rng;
};

class constant_dist: public rng_distribution {
  public:
	explicit constant_dist(float value)
		: m_value(value) { }

	float operator()()
	{
		return m_value;
	}
  private:
	float m_value;
};

class uniform_dist: public rng_distribution {
  public:
	explicit uniform_dist(float a, float b)
	{
		distribution = std::uniform_real_distribution<float>(a, b);
	}

	float operator()()
	{
		float value = 0.0;
		while (value < min_size)
			value = distribution(rng);
		return value;
	}

  private:
	std::uniform_real_distribution<float> distribution;
};

class normal_dist: public rng_distribution {
  public:
	explicit normal_dist(float mean, float sigma)
	{
		distribution = std::normal_distribution<float>(mean, sigma);
	}

	float operator()()
	{
		float value = 0.0;
		while (value < min_size)
			value = distribution(rng);
		return value;
	}

  private:
	std::normal_distribution<float> distribution;
};

class lognormal_dist: public rng_distribution {
  public:
	explicit lognormal_dist(float mean, float sigma)
	{
		distribution = std::lognormal_distribution<float>(mean, sigma);
	}

	float operator()()
	{
		float value = 0.0;
		while (value < min_size)
			value = distribution(rng);
		return value;
	}

  private:
	std::lognormal_distribution<float> distribution;
};

#if 0
class gamma_dist: public rng_distribution {
  public:
	explicit gamma_dist(float alpha, float beta)
	{
		distribution = std::gamma_distribution<float>(alpha, beta);
	}

	float operator()()
	{
		float value = 0.0;
		while (value < min_size)
			value = distribution(rng);
		return value;
	}

  private:
	std::gamma_distribution<float> distribution;
};
#endif

class cauchy_dist: public rng_distribution {
  public:
	explicit cauchy_dist(float a, float b)
	{
		distribution = std::cauchy_distribution<float>(a, b);
	}

	float operator()()
	{
		float value = 0.0;
		while (value < min_size)
			value = distribution(rng);
		return value;
	}

  private:
	std::cauchy_distribution<float> distribution;
};

class weibull_dist: public rng_distribution {
  public:
	explicit weibull_dist(float a, float b)
	{
		distribution = std::weibull_distribution<float>(a, b);
	}

	float operator()()
	{
		float value = 0.0;
		while (value < min_size)
			value = distribution(rng);
		return value;
	}

  private:
	std::weibull_distribution<float> distribution;
};

#endif
