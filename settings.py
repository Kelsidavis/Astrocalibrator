import os
import json

# Define full settings path inside AppData
SETTINGS_FILE = os.path.join(os.getenv('APPDATA'), 'Astrocalibrator', 'settings.json')
os.makedirs(os.path.dirname(SETTINGS_FILE), exist_ok=True)

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
