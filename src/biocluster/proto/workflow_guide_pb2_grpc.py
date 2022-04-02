# -*- coding: utf-8 -*-
# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
import grpc

import public_pb2 as public__pb2
import workflow_guide_pb2 as workflow__guide__pb2


class WorkflowGuideStub(object):
  """定义Workflow运行过程中与WFM的通信协议 
  """

  def __init__(self, channel):
    """Constructor.

    Args:
      channel: A grpc.Channel.
    """
    self.GetToolsStates = channel.unary_stream(
        '/workflowguid.WorkflowGuide/GetToolsStates',
        request_serializer=public__pb2.Workflow.SerializeToString,
        response_deserializer=workflow__guide__pb2.State.FromString,
        )
    self.StopTimeoutCheck = channel.unary_unary(
        '/workflowguid.WorkflowGuide/StopTimeoutCheck',
        request_serializer=public__pb2.Workflow.SerializeToString,
        response_deserializer=public__pb2.Success.FromString,
        )
    self.SendStep = channel.unary_unary(
        '/workflowguid.WorkflowGuide/SendStep',
        request_serializer=workflow__guide__pb2.Step.SerializeToString,
        response_deserializer=public__pb2.Success.FromString,
        )
    self.Update = channel.unary_unary(
        '/workflowguid.WorkflowGuide/Update',
        request_serializer=workflow__guide__pb2.Status.SerializeToString,
        response_deserializer=public__pb2.Success.FromString,
        )
    self.GetRunInfo = channel.unary_unary(
        '/workflowguid.WorkflowGuide/GetRunInfo',
        request_serializer=public__pb2.Workflow.SerializeToString,
        response_deserializer=workflow__guide__pb2.RunConfig.FromString,
        )
    self.AddBatch = channel.unary_unary(
        '/workflowguid.WorkflowGuide/AddBatch',
        request_serializer=workflow__guide__pb2.Batch.SerializeToString,
        response_deserializer=workflow__guide__pb2.BatchID.FromString,
        )
    self.Submit = channel.unary_unary(
        '/workflowguid.WorkflowGuide/Submit',
        request_serializer=workflow__guide__pb2.Job.SerializeToString,
        response_deserializer=workflow__guide__pb2.JobSuccess.FromString,
        )
    self.Delete = channel.unary_unary(
        '/workflowguid.WorkflowGuide/Delete',
        request_serializer=workflow__guide__pb2.Tool.SerializeToString,
        response_deserializer=public__pb2.Success.FromString,
        )


class WorkflowGuideServicer(object):
  """定义Workflow运行过程中与WFM的通信协议 
  """

  def GetToolsStates(self, request, context):
    """获取最新的与workflow相关的tool状态信息,当收到WFM的signal信号时执行
    当WFM多次通知，但是workflow无响应时，可以任务workflow被卡死
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def StopTimeoutCheck(self, request, context):
    """WFM对于超过2小时没有任何通信发生,且没有Tool或子任务运行的workflow(非即时任务),认为其发生死循环等异常情况,将强制终止
    对于进行IO阻塞需要超过2小时的情况(如大规模导表等),则需要在workflow中执行stop_timeout_check方法来阻止WFM强制终止进程情况的发生
    StopTimeout告诉WFM,停止对该workflow进行超时检测和强制终止
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def SendStep(self, request, context):
    """发送步骤信息
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def Update(self, request, context):
    """更新状态
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def GetRunInfo(self, request, context):
    """获取运行配置
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def AddBatch(self, request, context):
    """提交Batch并返回BatchID,对于重运行的任务,自动根据Path和参数匹配到原ID,并返回原BatchID
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def Submit(self, request, context):
    """提交tool任务    
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def Delete(self, request, context):
    """删除tool任务
    """
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')


def add_WorkflowGuideServicer_to_server(servicer, server):
  rpc_method_handlers = {
      'GetToolsStates': grpc.unary_stream_rpc_method_handler(
          servicer.GetToolsStates,
          request_deserializer=public__pb2.Workflow.FromString,
          response_serializer=workflow__guide__pb2.State.SerializeToString,
      ),
      'StopTimeoutCheck': grpc.unary_unary_rpc_method_handler(
          servicer.StopTimeoutCheck,
          request_deserializer=public__pb2.Workflow.FromString,
          response_serializer=public__pb2.Success.SerializeToString,
      ),
      'SendStep': grpc.unary_unary_rpc_method_handler(
          servicer.SendStep,
          request_deserializer=workflow__guide__pb2.Step.FromString,
          response_serializer=public__pb2.Success.SerializeToString,
      ),
      'Update': grpc.unary_unary_rpc_method_handler(
          servicer.Update,
          request_deserializer=workflow__guide__pb2.Status.FromString,
          response_serializer=public__pb2.Success.SerializeToString,
      ),
      'GetRunInfo': grpc.unary_unary_rpc_method_handler(
          servicer.GetRunInfo,
          request_deserializer=public__pb2.Workflow.FromString,
          response_serializer=workflow__guide__pb2.RunConfig.SerializeToString,
      ),
      'AddBatch': grpc.unary_unary_rpc_method_handler(
          servicer.AddBatch,
          request_deserializer=workflow__guide__pb2.Batch.FromString,
          response_serializer=workflow__guide__pb2.BatchID.SerializeToString,
      ),
      'Submit': grpc.unary_unary_rpc_method_handler(
          servicer.Submit,
          request_deserializer=workflow__guide__pb2.Job.FromString,
          response_serializer=workflow__guide__pb2.JobSuccess.SerializeToString,
      ),
      'Delete': grpc.unary_unary_rpc_method_handler(
          servicer.Delete,
          request_deserializer=workflow__guide__pb2.Tool.FromString,
          response_serializer=public__pb2.Success.SerializeToString,
      ),
  }
  generic_handler = grpc.method_handlers_generic_handler(
      'workflowguid.WorkflowGuide', rpc_method_handlers)
  server.add_generic_rpc_handlers((generic_handler,))


class NtmJobStub(object):
  """给NTM发送即时tool任务
  """

  def __init__(self, channel):
    """Constructor.

    Args:
      channel: A grpc.Channel.
    """
    self.Submit = channel.unary_unary(
        '/workflowguid.NtmJob/Submit',
        request_serializer=workflow__guide__pb2.Job.SerializeToString,
        response_deserializer=workflow__guide__pb2.JobSuccess.FromString,
        )
    self.StopWork = channel.unary_unary(
        '/workflowguid.NtmJob/StopWork',
        request_serializer=public__pb2.Workflow.SerializeToString,
        response_deserializer=public__pb2.Success.FromString,
        )


class NtmJobServicer(object):
  """给NTM发送即时tool任务
  """

  def Submit(self, request, context):
    # missing associated documentation comment in .proto file
    pass
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')

  def StopWork(self, request, context):
    # missing associated documentation comment in .proto file
    pass
    context.set_code(grpc.StatusCode.UNIMPLEMENTED)
    context.set_details('Method not implemented!')
    raise NotImplementedError('Method not implemented!')


def add_NtmJobServicer_to_server(servicer, server):
  rpc_method_handlers = {
      'Submit': grpc.unary_unary_rpc_method_handler(
          servicer.Submit,
          request_deserializer=workflow__guide__pb2.Job.FromString,
          response_serializer=workflow__guide__pb2.JobSuccess.SerializeToString,
      ),
      'StopWork': grpc.unary_unary_rpc_method_handler(
          servicer.StopWork,
          request_deserializer=public__pb2.Workflow.FromString,
          response_serializer=public__pb2.Success.SerializeToString,
      ),
  }
  generic_handler = grpc.method_handlers_generic_handler(
      'workflowguid.NtmJob', rpc_method_handlers)
  server.add_generic_rpc_handlers((generic_handler,))