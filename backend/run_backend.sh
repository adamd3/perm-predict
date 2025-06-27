source venv/bin/activate
uvicorn app.main:app --reload --port 8000

# Access the API at http://localhost:8000
# Interactive API docs: http://localhost:8000/docs

# Health check
curl http://localhost:8000/health
# {"status":"healthy","model_loaded":true,"model_version":"1.0.0","model_loaded_at":"2025-01-31T18:24:12.590458","timestamp":"2025-01-31T18:24:12.590467"}Adams-MacBook-Pro-8:backend adamdinan$ 
# Single prediction
curl -X POST http://localhost:8000/predict/single \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}'
# Internal Server Error

# Error log
# ERROR:    Exception in ASGI application
#   + Exception Group Traceback (most recent call last):
#   |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_utils.py", line 76, in collapse_excgroups
#   |     yield
#   |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 186, in __call__
#   |     async with anyio.create_task_group() as task_group:
#   |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/anyio/_backends/_asyncio.py", line 767, in __aexit__
#   |     raise BaseExceptionGroup(
#   | ExceptionGroup: unhandled errors in a TaskGroup (1 sub-exception)
#   +-+---------------- 1 ----------------
#     | Traceback (most recent call last):
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/uvicorn/protocols/http/h11_impl.py", line 403, in run_asgi
#     |     result = await app(  # type: ignore[func-returns-value]
#     |              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/uvicorn/middleware/proxy_headers.py", line 60, in __call__
#     |     return await self.app(scope, receive, send)
#     |            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/fastapi/applications.py", line 1054, in __call__
#     |     await super().__call__(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/applications.py", line 113, in __call__
#     |     await self.middleware_stack(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/errors.py", line 187, in __call__
#     |     raise exc
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/errors.py", line 165, in __call__
#     |     await self.app(scope, receive, _send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 185, in __call__
#     |     with collapse_excgroups():
#     |   File "/opt/local/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/contextlib.py", line 158, in __exit__
#     |     self.gen.throw(typ, value, traceback)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_utils.py", line 82, in collapse_excgroups
#     |     raise exc
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 187, in __call__
#     |     response = await self.dispatch_func(request, call_next)
#     |                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/app/main.py", line 26, in dispatch
#     |     response = await call_next(request)
#     |                ^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 163, in call_next
#     |     raise app_exc
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 149, in coro
#     |     await self.app(scope, receive_or_disconnect, send_no_error)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/exceptions.py", line 62, in __call__
#     |     await wrap_app_handling_exceptions(self.app, conn)(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 53, in wrapped_app
#     |     raise exc
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 42, in wrapped_app
#     |     await app(scope, receive, sender)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 715, in __call__
#     |     await self.middleware_stack(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 735, in app
#     |     await route.handle(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 288, in handle
#     |     await self.app(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 76, in app
#     |     await wrap_app_handling_exceptions(app, request)(scope, receive, send)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 53, in wrapped_app
#     |     raise exc
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 42, in wrapped_app
#     |     await app(scope, receive, sender)
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 73, in app
#     |     response = await f(request)
#     |                ^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/fastapi/routing.py", line 301, in app
#     |     raw_response = await run_endpoint_function(
#     |                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/fastapi/routing.py", line 212, in run_endpoint_function
#     |     return await dependant.call(**values)
#     |            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/slowapi/extension.py", line 734, in async_wrapper
#     |     response = await func(*args, **kwargs)  # type: ignore
#     |                ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/app/main.py", line 59, in predict_single
#     |     results = await process_batch([input_data.smiles], session)
#     |               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/app/utils/processing.py", line 38, in process_batch
#     |     predictions = session.run(None, {input_name: features_array})[0]
#     |                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     |   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/onnxruntime/capi/onnxruntime_inference_collection.py", line 266, in run
#     |     return self._sess.run(output_names, input_feed, run_options)
#     |            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     | onnxruntime.capi.onnxruntime_pybind11_state.InvalidArgument: [ONNXRuntimeError] : 2 : INVALID_ARGUMENT : Unexpected input data type. Actual: (tensor(int8)) , expected: (tensor(float))
#     +------------------------------------

# During handling of the above exception, another exception occurred:

# Traceback (most recent call last):
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/uvicorn/protocols/http/h11_impl.py", line 403, in run_asgi
#     result = await app(  # type: ignore[func-returns-value]
#              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/uvicorn/middleware/proxy_headers.py", line 60, in __call__
#     return await self.app(scope, receive, send)
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/fastapi/applications.py", line 1054, in __call__
#     await super().__call__(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/applications.py", line 113, in __call__
#     await self.middleware_stack(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/errors.py", line 187, in __call__
#     raise exc
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/errors.py", line 165, in __call__
#     await self.app(scope, receive, _send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 185, in __call__
#     with collapse_excgroups():
#   File "/opt/local/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/contextlib.py", line 158, in __exit__
#     self.gen.throw(typ, value, traceback)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_utils.py", line 82, in collapse_excgroups
#     raise exc
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 187, in __call__
#     response = await self.dispatch_func(request, call_next)
#                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/app/main.py", line 26, in dispatch
#     response = await call_next(request)
#                ^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 163, in call_next
#     raise app_exc
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/base.py", line 149, in coro
#     await self.app(scope, receive_or_disconnect, send_no_error)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/middleware/exceptions.py", line 62, in __call__
#     await wrap_app_handling_exceptions(self.app, conn)(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 53, in wrapped_app
#     raise exc
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 42, in wrapped_app
#     await app(scope, receive, sender)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 715, in __call__
#     await self.middleware_stack(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 735, in app
#     await route.handle(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 288, in handle
#     await self.app(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 76, in app
#     await wrap_app_handling_exceptions(app, request)(scope, receive, send)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 53, in wrapped_app
#     raise exc
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/_exception_handler.py", line 42, in wrapped_app
#     await app(scope, receive, sender)
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/starlette/routing.py", line 73, in app
#     response = await f(request)
#                ^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/fastapi/routing.py", line 301, in app
#     raw_response = await run_endpoint_function(
#                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/fastapi/routing.py", line 212, in run_endpoint_function
#     return await dependant.call(**values)
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/slowapi/extension.py", line 734, in async_wrapper
#     response = await func(*args, **kwargs)  # type: ignore
#                ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/app/main.py", line 59, in predict_single
#     results = await process_batch([input_data.smiles], session)
#               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/app/utils/processing.py", line 38, in process_batch
#     predictions = session.run(None, {input_name: features_array})[0]
#                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/Users/adamdinan/perm-predict/backend/venv/lib/python3.11/site-packages/onnxruntime/capi/onnxruntime_inference_collection.py", line 266, in run
#     return self._sess.run(output_names, input_feed, run_options)
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# onnxruntime.capi.onnxruntime_pybind11_state.InvalidArgument: [ONNXRuntimeError] : 2 : INVALID_ARGUMENT : Unexpected input data type. Actual: (tensor(int8)) , expected: (tensor(float))

# Batch prediction (with a CSV file)
curl -X POST http://localhost:8000/predict/batch \
  -F "file=@your_smiles.csv"