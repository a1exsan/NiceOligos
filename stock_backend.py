from OligoMap_utils import api_db_interface

class stock_backend_model(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)
        self.pincode = ''