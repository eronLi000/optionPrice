import re

def scan_file(path_to_file):
    all_data = {}
    with open(path_to_file, 'r') as eFile:
        for line in eFile:
            _stockInfo = line.split(',', 1)
            if re.search(r'^'+""+'(.*?)$', _stockInfo[0], re.IGNORECASE):
                all_data[_stockInfo[0]] = _stockInfo[1].strip('\n')
    eFile.close()
    return all_data