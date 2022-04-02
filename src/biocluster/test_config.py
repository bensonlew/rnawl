# -*- coding: utf-8 -*-

import os
import socket
import grpc
from .proto import public_pb2, config_pb2, config_pb2_grpc

def main():
    port = 7321
    hostname = socket.gethostname()
    ip = socket.gethostbyname(hostname)
    # print("获取webauth:hostname %s, ip %s,port %s",hostname,ip,port)
    print "获取webauth:hostname {}, ip {},port {}".format(hostname,ip,port)
    try:
        with grpc.insecure_channel('localhost:%s' % port) as channel:
            stub = config_pb2_grpc.ConfigServerStub(channel)
            response = stub.GetWebAuth(config_pb2.ProjectType(
                type="tool_lab",
            ))
            # authkey: 458b97de4c0bb5bf416c8cea208309ed
            # client: client03
            # binds_id: 5e8c0a091b1800007c006b1a
            # interface_id: 1348
            # env_name: offline
            cfg = {
                "authkey": response.key,
                "client": response.client,
                "binds_id": response.bindsid,
                "env_name": response.envname,
            }
            self._webauth_info[ptype] = cfg
            self._grpc_times = 0
    except Exception as e:
        exstr = traceback.format_exc()
        print(exstr)
        sys.stdout.flush()

if __name__ == '__main__':
    main()