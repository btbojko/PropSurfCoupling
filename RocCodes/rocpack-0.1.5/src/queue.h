#ifndef _QUEUE_H
#define _QUEUE_H

#include <vector>

#include "event.h"

class EventQueue {
  public:
	EventQueue(unsigned int max_size);

	Event* next_event() const;

	void push(Event* event);
	void pop();
	void make_heap();
	void update();
	void update_id(unsigned int id);
	void update_all();

  private:
	std::vector<Event*> heap;
	std::vector<unsigned int>index;

	void swap(unsigned int i, unsigned int j);
	void move_up(unsigned int n);
	void move_down(unsigned int n);
	void update_event(Event* event);
};

extern EventQueue* event_queue;

#endif
