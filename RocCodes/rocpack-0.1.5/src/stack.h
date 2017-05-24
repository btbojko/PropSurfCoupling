#ifndef _STACK_H
#define _STACK_H

#include "scalar.h"

class Particle;

class Stack {
  public:
	Stack(int max_elements) {
		n = 0; max = max_elements;
		m_stack = (Particle **) malloc(max_elements * sizeof(Particle*));
	}

	~Stack() { free(m_stack); }

	void push(Particle* p)
	{
		assert(n < max);
		m_stack[n++] = p;
	}

	Particle* top()
	{
		assert(n > 0);
		return m_stack[n-1];
	}

	void pop()
	{
		assert(n > 0);
		n--;
	}

	void clear()
	{
		n = 0;
	}

	int size()
	{
		return n;
	}

	bool empty()
	{
		return n == 0;
	}

  private:
	int n, max;
	Particle** m_stack;
};

#endif

