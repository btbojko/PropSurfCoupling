#include <vector>

#include "event.h"
#include "queue.h"

EventQueue* event_queue;

EventQueue::EventQueue(unsigned int max_size)
{
	heap.reserve(max_size);
	index.resize(max_size);
}

Event* EventQueue::next_event() const { return heap.front(); }

void EventQueue::push(Event* event)
{
	heap.push_back(event);
	index[event->id()] = heap.size() - 1;
	move_up(heap.size() - 1);
}

void EventQueue::pop()
{
	/* swap front and back */
	heap.front() = heap.back();
	index[heap.back()->id() ] = 0;
	heap.pop_back();
	move_down(0);
}

void EventQueue::make_heap()
{
	for (unsigned int i = 0; i < heap.size(); i++)
		move_up(i);
}

void EventQueue::update()
{
	if (heap.front()->type() == COLLISION) {
		unsigned int id = heap.front()->secondary_id();
		heap.front()->recompute();
		move_down(0);
		heap[ index[id] ]->recompute();
		move_up(index[id]);
		move_down(index[id]);
	} else {
		heap.front()->recompute();
		move_down(0);
	}
}

void EventQueue::update_id(unsigned int id)
{
	unsigned int pid = index[id];
	if (heap[pid]->type() == COLLISION) {
		unsigned int sid = heap[pid]->secondary_id();
		heap[pid]->recompute();
		heap[sid]->recompute();
		move_up(pid);
		move_down(pid);
		move_up(sid);
		move_down(sid);
	} else {
		heap[pid]->recompute();
		move_up(pid);
		move_down(pid);
	}
}

void EventQueue::update_all()
{
	unsigned int i;
	for (i = 0; i < heap.size(); i++) {
		heap[i]->recompute();
	}

	make_heap();
}

void EventQueue::swap(unsigned int i, unsigned int j)
{
	Event* event;
	index[ heap[i]->id() ] = j;
	index[ heap[j]->id() ] = i;
	event   = heap[i];
	heap[i] = heap[j];
	heap[j] = event;
}

void EventQueue::move_up(unsigned int n)
{
	if (n == 0) return;

	unsigned int p;

	do {
		p = (n-1)/2;
		if (heap[n]->time() < heap[p]->time()) {
			swap(n, p); n = p;
		} else {
			break;
		}
	} while (p > 0);
}

void EventQueue::move_down(unsigned int n)
{
	unsigned int c = 2*n+1;
	while (c < heap.size()) {
		if (c + 1 < heap.size())
			if (heap[c+1]->time() < heap[c]->time())
				c = c + 1;

		if (heap[n]->time() > heap[c]->time()) {
			swap(n, c); n = c; c = 2*c+1;
		} else {
			break;
		}
	}
}

void EventQueue::update_event(Event* event)
{
	move_up(index[event->id()]);
	move_down(index[event->id()]);
}

