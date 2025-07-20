
from OligoMap_utils import api_db_interface


class Oligomap_backend(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        self.client = {}
        self.client_frontend = {}

    def init_frontend(self, front, ip):
        if ip not in list(self.client.keys()):
            self.client[ip] = ''
            self.client_frontend[ip] = front.get_model()
        self.frontend = front
        self.frontend.set_model(self.client_frontend[ip])
        #print(self.client_frontend[ip])