#ifndef _SETTINGS_H
#define _SETTINGS_H

namespace settings {
	bool isset(const char *name);
	template <typename T> T get(const char *name);
	bool isset(const char *name);
	void set(const char *name, bool value);
	void set(const char *name, int value);
	void set(const char *name, float value);
	void set(const char *name, double value);
	void set(const char *name, const char *value);
	void unset(const char *name);
}

#endif
