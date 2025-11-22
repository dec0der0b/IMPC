import requests
REQUESTS_URL = 'https://rest.xialab.ca/api/mapcompounds'


class MetaboAnalyst:

    @staticmethod
    def query(metabolites: list[str]):
        payload = '{"queryList": "{' + ';'.join(metabolites) + '"}","inputType": "name"}'
        headers = {
            'Content-Type': "application/json",
            'cache-control': "no-cache",
        }

        response = requests.request("POST", REQUESTS_URL, data=payload, headers=headers)

        return response.text

