# pytest fixtures for ascairn tests
import os
import subprocess
import pytest

# Limit Polars threads to avoid memory issues during tests
os.environ["POLARS_MAX_THREADS"] = "1"

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="session")
def resource_dir():
    """Clone ascairn_resource (devel branch) if not exists and return path."""
    resource_path = os.path.join(TESTS_DIR, "ascairn_resource")
    if not os.path.exists(resource_path):
        subprocess.run(
            ["git", "clone", "-b", "devel", "https://github.com/friend1ws/ascairn_resource.git"],
            cwd=TESTS_DIR,
            check=True,
        )
    return os.path.join(resource_path, "resource", "ver_2024-12-06")


@pytest.fixture(scope="session")
def output_dir():
    path = os.path.join(TESTS_DIR, "output")
    os.makedirs(path, exist_ok=True)
    return path


@pytest.fixture
def expected_dir():
    return os.path.join(TESTS_DIR, "expected")


@pytest.fixture
def s3_cram():
    return "s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram"
