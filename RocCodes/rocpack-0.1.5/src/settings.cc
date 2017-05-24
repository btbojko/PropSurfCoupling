#include <map>
#include <string>

namespace settings {
	static std::map<std::string, std::string> settings;

	/* get methods */

	template <typename T> T get(const char *name);

	bool isset(const char *name)
	{
		return settings.find(name) != settings.end();
	}

	template <> bool get(const char *name)
	{
		return settings[name] == "true";
	}

	template <> int   get(const char *name)
	{
		return stoi(settings[name]);
	}

	template <> float get(const char *name)
	{
		return stof(settings[name]);
	}

	template <> double get(const char *name)
	{
		return stod(settings[name]);
	}

	template <> const char* get(const char *name)
	{
		return settings[name].c_str();
	}

	/* set methods */

	void set(const char *name, bool value)
	{
		settings[name] = (value ? "true" : "false");
	}

	void set(const char *name, int value)
	{
		settings[name] = std::to_string(value);
	}

	void set(const char *name, float value)
	{
		settings[name] = std::to_string(value);
	}

	void set(const char *name, double value)
	{
		settings[name] = std::to_string(value);
	}

	void set(const char *name, const char *value)
	{
		settings[name] = value;
	}

	void unset(const char *name)
	{
		settings.erase(name);
	}
}
