import subprocess
import platform


def check_tool(name: str, version_arg: str = "--version", timeout: float = 5.0):
    try:
        completed = subprocess.run(
            [name, version_arg],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            check=False,
        )
        return True, completed.stdout.strip()
    except FileNotFoundError:
        return False, None
    except Exception as e:
        return False, str(e)


def ensure_tool_installed(name: str):
    installed, _ = check_tool(name)
    if not installed:
        raise OSError(f"{name} is not installed.")
    return True


def ensure_unix():
    if platform.system() == "Windows":
        raise OSError("This operation is not supported on Windows.")
    else:
        return True
