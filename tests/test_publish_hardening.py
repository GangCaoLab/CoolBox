import subprocess as subp
import sys
import zipfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PUBLISH_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "python-publish.yml"
DIST_CHECKER = REPO_ROOT / "scripts" / "check_dist_safety.py"


def run_dist_checker(*paths):
    return subp.run(
        [sys.executable, str(DIST_CHECKER), *map(str, paths)],
        cwd=REPO_ROOT,
        text=True,
        stdout=subp.PIPE,
        stderr=subp.STDOUT,
    )


def make_wheel(path, entries):
    with zipfile.ZipFile(path, "w") as wheel:
        for name, content in entries.items():
            wheel.writestr(name, content)


def test_python_publish_uses_trusted_publishing_without_long_lived_tokens():
    workflow = PUBLISH_WORKFLOW.read_text()

    assert "pypa/gh-action-pypi-publish@release/v1" in workflow
    assert "id-token: write" in workflow
    assert "secrets.PYPI_API_TOKEN" not in workflow
    assert "twine upload" not in workflow


def test_dist_checker_accepts_clean_wheel(tmp_path):
    wheel = tmp_path / "coolbox-0.4.3-py3-none-any.whl"
    make_wheel(
        wheel,
        {
            "coolbox/__init__.py": "__version__ = '0.4.3'\n",
            "coolbox-0.4.3.dist-info/METADATA": "Name: coolbox\n",
            "coolbox-0.4.3.dist-info/RECORD": "",
        },
    )

    result = run_dist_checker(wheel)

    assert result.returncode == 0, result.stdout


def test_dist_checker_allows_normal_python_subprocess_usage(tmp_path):
    wheel = tmp_path / "coolbox-0.4.3-py3-none-any.whl"
    make_wheel(
        wheel,
        {
            "coolbox/utilities/cmd.py": "import subprocess\nsubprocess.run(['true'])\n",
            "coolbox-0.4.3.dist-info/METADATA": "Name: coolbox\n",
            "coolbox-0.4.3.dist-info/RECORD": "",
        },
    )

    result = run_dist_checker(wheel)

    assert result.returncode == 0, result.stdout


def test_dist_checker_rejects_executable_pth_and_javascript_payload(tmp_path):
    wheel = tmp_path / "coolbox-0.4.3-py3-none-any.whl"
    make_wheel(
        wheel,
        {
            "coolbox/__init__.py": "__version__ = '0.4.3'\n",
            "coolbox-setup.pth": "import os\n",
            "coolbox/_index.js": "console.log('unexpected payload')\n",
            "coolbox-0.4.3.dist-info/METADATA": "Name: coolbox\n",
            "coolbox-0.4.3.dist-info/RECORD": "",
        },
    )

    result = run_dist_checker(wheel)

    assert result.returncode != 0
    assert "coolbox-setup.pth" in result.stdout
    assert "coolbox/_index.js" in result.stdout
