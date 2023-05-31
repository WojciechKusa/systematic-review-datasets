import argparse
import json
import os

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config/")
CONFIG_FILES = {
    "s2_api_key": "semanticscholar_api_key.txt",
    "core_api_key": "core_api_key.txt",
    "entrez_email": "entrez_email.txt",
    "cochrane_cookie": "cochrane_cookie.txt",
    "grobid_config": "grobid_config.json",
}

grobid_config = {
    "batch_size": 1000,
    "sleep_time": 5,
    "timeout": 60,
    "coordinates": ["persName", "figure", "ref", "biblStruct", "formula", "s"],
}


def _update_grobid_config() -> None:
    grobid_server = input(
        "Enter your Grobid server (e.g., 'http://192.192.111.12:8070'): "
    )
    if grobid_server:
        grobid_config["grobid_server"] = grobid_server
        with open(os.path.join(CONFIG_PATH, CONFIG_FILES["grobid_config"]), "w+") as f:
            json.dump(grobid_config, f, indent=4)


def configure(overwrite: bool, update: str) -> None:
    """
    Configure the environment variables. This is necessary for the scripts to run.
    """
    if not os.path.exists(CONFIG_PATH):
        os.mkdir(CONFIG_PATH)

    for config_key, config_file in CONFIG_FILES.items():
        file_path = os.path.join(CONFIG_PATH, config_file)
        if not overwrite and os.path.exists(file_path) and update != config_key:
            print(
                f"{config_key} already exists. Use --overwrite or --update to change."
            )
            continue

        if config_key == "grobid_config":
            _update_grobid_config()
        else:
            value = input(f"Enter your {config_key.replace('_', '.')} value: ")
            if value:
                with open(file_path, "w+") as f:
                    f.write(value)

    print("Configuration complete. You can now run the scripts.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite all existing configurations"
    )
    parser.add_argument(
        "--update",
        type=str,
        choices=CONFIG_FILES.keys(),
        help="Update a specific configuration",
        default=None,
    )
    args = parser.parse_args()

    configure(overwrite=args.overwrite, update=args.update)
