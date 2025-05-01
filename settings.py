import os
import sys
import json

# Cross-platform settings directory function
def get_app_data_dir():
    # Windows: use APPDATA
    if os.name == 'nt' and os.getenv('APPDATA'):
        return os.path.join(os.getenv('APPDATA'), 'Astrocalibrator')
    # macOS: use ~/Library/Application Support
    elif sys.platform == 'darwin':
        return os.path.join(os.path.expanduser('~'), 'Library', 'Application Support', 'Astrocalibrator')
    # Linux/Unix: use XDG_CONFIG_HOME or ~/.config
    else:
        config_home = os.getenv('XDG_CONFIG_HOME') or os.path.join(os.path.expanduser('~'), '.config')
        return os.path.join(config_home, 'Astrocalibrator')

# Create settings directory and file path
SETTINGS_DIR = get_app_data_dir()
SETTINGS_FILE = os.path.join(SETTINGS_DIR, 'settings.json')
os.makedirs(SETTINGS_DIR, exist_ok=True)

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
