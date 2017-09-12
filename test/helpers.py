import os.path

def get_test_data_dir():
    return os.path.join(os.path.dirname(__file__), '..', 'test_data')
