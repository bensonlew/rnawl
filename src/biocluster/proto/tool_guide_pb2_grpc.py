# -*- coding: utf-8 -*-
# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
import grpc

import tool_guide_pb2 as tool__guide__pb2


class ToolGuideStub(object):
  """tool运行期间与ntm的通信
  """

  def __init__(self, channel):
    """Constructor.

    Args:
      channel: A grpc.Channel.
    """
    self.SendState = channel.stream_unary(
        '/ntmpb.ToolGuide/SendState',
        request_serializer=tool__guide__pb2.State.SerializeToString,
        response_deserializer=tool__guide__pb2.Success.FromString,
        )
    self.OptionData = channel.unary_unary(
        '/ntmpb.ToolGuide/OptionData',
        request_serializer=tool__guide__pb2.Options.SerializeToString,
        response_deserializer=tool__guide__pb2.Success.FromString,
        )
    self.Keepalive = channel.unary_unary(
        '/ntmpb.ToolGuide/Keepalive',
        request_serializer=tool__guide__pb2.Tool.SerializeToString,
        response_deserializer=tool__guide__pb2.Success.FromString,
        )


class ToolGuideServicer(object):
  """tool运行期间与ntm的通信
  """

  def SendState(self, request_iterator, context):
    """发送状态
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def OptionData(self, request, context):
    """发送参数信息
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def Keepalive(self, request, context):
    """keepalive
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')


def add_ToolGuideServicer_to_server(servicer, server):
  rpc_method_handlers = {
      'SendState': grpc.stream_unary_rpc_method_handler(
          servicer.SendState,
          request_deserializer=tool__guide__pb2.State.FromString,
          response_serializer=tool__guide__pb2.Success.SerializeToString,
      ),
      'OptionData': grpc.unary_unary_rpc_method_handler(
          servicer.OptionData,
          request_deserializer=tool__guide__pb2.Options.FromString,
          response_serializer=tool__guide__pb2.Success.SerializeToString,
      ),
      'Keepalive': grpc.unary_unary_rpc_method_handler(
          servicer.Keepalive,
          request_deserializer=tool__guide__pb2.Tool.FromString,
          response_serializer=tool__guide__pb2.Success.SerializeToString,
      ),
  }
  generic_handler = grpc.method_handlers_generic_handler(
      'ntmpb.ToolGuide', rpc_method_handlers)
  server.add_generic_rpc_handlers((generic_handler,))
