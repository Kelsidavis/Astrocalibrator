import os
import json

SETTINGS_FILE = os.path.join(os.path.expanduser('~'), '.astro_calibrator_settings.json')

def load_settings():
    try:
        with open(SETTINGS_FILE, 'r') as f:
            return json.load(f)
    except Exception:
        return {}

def save_settings(data):
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(data, f, indent=4)
    except Exception as e:
        print(f"⚠️ Failed to save settings: {e}")

def get_remembered_file(key):
    settings = load_settings()
    return settings.get(key)

def remember_file(key, path):
    settings = load_settings()
    settings[key] = path
    save_settings(settings)
