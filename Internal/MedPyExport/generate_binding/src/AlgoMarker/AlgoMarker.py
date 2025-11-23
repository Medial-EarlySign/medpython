import ctypes, json, traceback, os
from functools import wraps
from typing import Any, List, Callable


class SingleDataElement:
    """SingleDataElement object that holds a single data element in AlgoMarker patient data repository."""

    times: List[int]
    values: List[float]

    def __init__(self, times: List[int], values: List[float]):
        """SingleDataElement constructor - receives signal name, times and values"""
        self.times = times
        self.values = values

    def __repr__(self):
        return f"(times={self.times}, values={self.values})"


class AlgoMarker:
    """AlgoMarker object that holds full model pipeline to calculate meaningfull insights from EMR raw data.

    Methods
    -------
    calculate
        recieves a request for execution of the model pipeline and returns a responde
    discovery
        returns a json specification of the AlgoMarker information, inputs, etc.
    dispose
        Release object memory - recomanded to use "with" statement
    clear_data
        clears AlgoMarker patient data repository memory
    add_data
        loads the AlgoMarker patient data repository memory with patient data
    """

    def __test_not_disposed(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args):
            s_obj = args[0]
            if s_obj.__disposed:
                raise NameError(
                    f"Error - Can't call {func.__name__} after algomarker was disposed"
                )
            return func(*args)

        return wrapper

    @staticmethod
    def __load_am_lib(libpath: str) -> tuple[ctypes.CDLL, int]:
        api_level = 2
        # Load the shared library into ctypes
        c_lib = ctypes.CDLL(libpath)
        c_lib.AM_API_Create.argtypes = (ctypes.c_int32, ctypes.POINTER(ctypes.c_void_p))
        c_lib.AM_API_Load.argtypes = (ctypes.c_void_p, ctypes.POINTER(ctypes.c_char))
        c_lib.AM_API_Load.restype = ctypes.c_int32
        c_lib.AM_API_DisposeAlgoMarker.argtypes = [ctypes.c_void_p]
        c_lib.AM_API_DisposeAlgoMarker.restype = None
        c_lib.AM_API_ClearData.argtypes = [ctypes.c_void_p]

        if (
            hasattr(c_lib, "AM_API_AddDataByType")
            and hasattr(c_lib, "AM_API_CalculateByType")
            and hasattr(c_lib, "AM_API_Discovery")
        ):
            c_lib.AM_API_Discovery.argtypes = (
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_char_p),
            )
            c_lib.AM_API_Discovery.restype = None
            c_lib.AM_API_AddDataByType.argtypes = (
                ctypes.c_void_p,
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_char_p),
            )
            c_lib.AM_API_CalculateByType.argtypes = (
                ctypes.c_void_p,
                ctypes.c_int32,
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_char_p),
            )
            c_lib.AM_API_Dispose.argtypes = [ctypes.c_char_p]
            c_lib.AM_API_Dispose.restype = None
        else:
            print(
                "Warning: AM_API_AddDataByType or AM_API_CalculateByType not found in the library, using old API"
            )
            api_level = 1

        c_lib.AM_API_AddData.argtypes = (
            ctypes.c_void_p,
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_char),
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_long),
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_float),
        )
        c_lib.AM_API_AddData.restype = ctypes.c_int32
        c_lib.AM_API_GetName.argtypes = (
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_char_p),
        )
        c_lib.AM_API_GetName.restype = None

        c_lib.AM_API_CreateRequest.argtypes = (
            ctypes.c_char_p,  # Request type
            ctypes.POINTER(ctypes.c_char_p),  # Score Type
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_long),
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_void_p),
        )
        c_lib.AM_API_CreateRequest.restype = ctypes.c_int32
        c_lib.AM_API_CreateResponses.argtypes = (
            ctypes.POINTER(ctypes.c_void_p),  # AlgoMarker object
        )  # Request type
        c_lib.AM_API_CreateResponses.restype = None
        c_lib.AM_API_DisposeRequest.argtypes = (ctypes.c_void_p,)  # Request object
        c_lib.AM_API_DisposeRequest.restype = None
        c_lib.AM_API_DisposeResponses.argtypes = (ctypes.c_void_p,)  # Response object
        c_lib.AM_API_DisposeResponses.restype = None
        c_lib.AM_API_Calculate.argtypes = (
            ctypes.c_void_p,  # AlgoMarker object
            ctypes.c_void_p,  # Request object
            ctypes.c_void_p,  # Response json string
        )
        c_lib.AM_API_Calculate.restype = ctypes.c_int32

        c_lib.AM_API_GetResponsesNum.argtypes = (ctypes.c_void_p,)
        c_lib.AM_API_GetResponsesNum.restype = ctypes.c_int32

        c_lib.AM_API_GetResponseAtIndex.argtypes = (
            ctypes.c_void_p,
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_void_p),
        )
        c_lib.AM_API_GetResponseAtIndex.restype = ctypes.c_int32

        c_lib.AM_API_GetResponseScoresNum.argtypes = (
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_int32),
        )
        c_lib.AM_API_GetResponseScoresNum.restype = ctypes.c_int32

        c_lib.AM_API_GetResponsePoint.argtypes = (
            ctypes.c_void_p,  # Response object
            ctypes.POINTER(ctypes.c_int32),  # Patient ID
            ctypes.POINTER(ctypes.c_long),  # Timestamp
        )
        c_lib.AM_API_GetResponsePoint.restype = ctypes.c_int32

        c_lib.AM_API_GetResponseMessages.argtypes = (
            ctypes.c_void_p,  # Response object
            ctypes.POINTER(ctypes.c_int32),  # Number of messages
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)),  # Message codes
            ctypes.POINTER(ctypes.POINTER(ctypes.c_char_p)),  # Messages errors
        )
        c_lib.AM_API_GetResponseMessages.restype = ctypes.c_int32

        c_lib.AM_API_GetScoreMessages.argtypes = (
            ctypes.c_void_p,  # Response object
            ctypes.c_int32,  # score_index
            ctypes.POINTER(ctypes.c_int32),  # Number of messages
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)),  # Message codes
            ctypes.POINTER(ctypes.POINTER(ctypes.c_char_p)),  # Messages errors
        )
        c_lib.AM_API_GetScoreMessages.restype = ctypes.c_int32

        c_lib.AM_API_GetResponseScoreByIndex.argtypes = (
            ctypes.c_void_p,  # Response object
            ctypes.c_int32,  # score_index
            ctypes.POINTER(ctypes.c_float),  # Score value
            ctypes.POINTER(ctypes.c_char_p),  # Score type
        )
        c_lib.AM_API_GetResponseScoreByIndex.restype = ctypes.c_int32

        c_lib.AM_API_GetSharedMessages.argtypes = (
            ctypes.c_void_p,  # Response object
            ctypes.POINTER(ctypes.c_int32),  # Number of messages
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)),  # Message codes
            ctypes.POINTER(ctypes.POINTER(ctypes.c_char_p)),  # Messages errors
        )
        c_lib.AM_API_GetSharedMessages.restype = ctypes.c_int32

        return c_lib, api_level

    @staticmethod
    def create_request_json(patient_id: int, prediction_time: int) -> str:
        """Creates and returns a string json request for patient_id and prediction_time"""
        js_req = (
            '{"type": "request", "request_id": "REQ_ID_1234", '
            + '"export": {"prediction": "pred_0"}, "requests": [ '
            + '{"patient_id":"%d", "time": "%d"} ]}'
            % (int(patient_id), int(prediction_time))
        )
        return js_req

    def __init__(self, amconfig_path: str, libpath: str | None = None):
        """AlgoMarker constractor - receives AlgoMarker configuration file path "amconfig".
        Optional path to C shared library file. If we want to use other version, not default
        library that is packed in this module.
        """
        if libpath is None:
            libpath = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "libdyn_AlgoMarker.so"
            )
        self.__lib = None
        self.__lib, self.api_version = AlgoMarker.__load_am_lib(libpath)
        self.__libpath = libpath
        print(f"Loaded library from {self.__libpath}")
        self.__obj = ctypes.c_void_p()
        res = self.__lib.AM_API_Create(1, ctypes.pointer(self.__obj))
        if res != 0:
            print("Error in creating AlgoMarker object")
        self.__disposed = False
        self.__name = None
        self.__amconfig_path = amconfig_path
        self.__load_algomarker(amconfig_path)

    def __load_algomarker(self, amconfig_path: str):
        if not (os.path.exists(amconfig_path)):
            raise NameError(
                f'amconfig path "{amconfig_path}" not found. File Not Found'
            )
        assert self.__lib is not None
        am_path = ctypes.create_string_buffer(amconfig_path.encode("ascii"))
        res = self.__lib.AM_API_Load(self.__obj, am_path)

        if res != 0:
            raise NameError(f"Error in loading AlgoMarker: {res}")
        else:
            try:
                info_js = self.discovery()
                if "name" in info_js:
                    self.__name = info_js["name"]
                    print(f"Loaded {self.__name} AlgoMarker succefully")
            except:
                print("Warning: couldn't retrieve AlgoMarker Name")

    def __repr__(self):
        if self.__disposed:
            return f"AlgoMarker was loaded with library {self.__libpath} and amconfig {self.__amconfig_path}, but disposed!"
        if self.__name is not None:
            return f"AlgoMarker {self.__name} was loaded with library {self.__libpath} and amconfig {self.__amconfig_path}"
        else:
            return f"AlgoMarker was loaded with library {self.__libpath} and amcofig {self.__amconfig_path}"

    def dispose(self):
        """Disposes the AlgoMarker object and frees the memory"""
        if self.__lib is not None:
            self.__lib.AM_API_DisposeAlgoMarker(self.__obj)
            self.__disposed = True
            if self.__name is None:
                print("Released AlgoMarker object")
            else:
                print(f'Released "{self.__name}" AlgoMarker object')
            self.__lib = None
            self.__obj = None

    def __del__(self):
        self.dispose()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.dispose()

    @__test_not_disposed
    def __dispose_string_mem(self, obj):
        assert self.__lib is not None
        self.__lib.AM_API_Dispose(obj)

    @__test_not_disposed
    def get_name(self) -> dict[str, Any]:
        """Returns information about the Algomarkers in json format - input signals, name, version, etc."""
        assert self.__lib is not None
        res_name = ctypes.c_char_p()
        self.__lib.AM_API_GetName(self.__obj, ctypes.byref(res_name))
        try:
            if res_name.value is None:
                raise NameError("Error in getting AlgoMarker name - name is None")
            res_discovery_str = res_name.value.decode("ascii")
            # Clear memory:
            res_discovery_str = {"name": res_discovery_str}
            return res_discovery_str
        except:
            print("Error in discovery json conversion")
            traceback.print_exc()
            raise

    @__test_not_disposed
    def discovery(self) -> dict[str, Any]:
        """Returns information about the Algomarkers in json format - input signals, name, version, etc."""
        assert self.__lib is not None
        if self.api_version == 1:
            return self.get_name()
        res_discovery = ctypes.c_char_p()
        self.__lib.AM_API_Discovery(self.__obj, ctypes.byref(res_discovery))
        try:
            res_discovery_str = res_discovery.value
            # Clear memory:
            self.__dispose_string_mem(res_discovery)
            if res_discovery_str is None:
                raise NameError(
                    "Error in getting AlgoMarker discovery - discovery is None"
                )
            res_discovery_str = json.loads(res_discovery_str)
            return res_discovery_str
        except:
            print("Error in discovery json conversion")
            traceback.print_exc()
            raise

    @__test_not_disposed
    def clear_data(self):
        """Frees the algomarker patient data repository"""
        assert self.__lib is not None
        res = self.__lib.AM_API_ClearData(self.__obj)
        if res != 0:
            raise NameError(f"Error in clearing data - error code {res}")

    @__test_not_disposed
    def add_data_simple(
        self, patient_id: int, signal_name: str, data: List[SingleDataElement]
    ) -> list[str]:
        assert self.__lib is not None
        flat_times = []
        flat_values = []
        times_size = None
        values_size = None
        messages = []
        for elem in data:
            flat_times.extend(elem.times)
            flat_values.extend(elem.values)
            if times_size is None:
                times_size = len(elem.times)
            if values_size is None:
                values_size = len(elem.values)
            if len(elem.times) != times_size:
                raise ValueError(
                    f"Error in add_data_simple - all times must have the same size, but got {len(elem.times)} != {times_size}"
                )
            if len(elem.values) != values_size:
                raise ValueError(
                    f"Error in add_data_simple - all values must have the same size, but got {len(elem.values)} != {values_size}"
                )

            # Convert to ctypes arrays
            c_times = (ctypes.c_long * len(flat_times))(*flat_times)
            c_values = (ctypes.c_float * len(flat_values))(*flat_values)
            res = self.__lib.AM_API_AddData(
                self.__obj,
                patient_id,
                ctypes.create_string_buffer(signal_name.encode("ascii")),
                len(flat_times),
                c_times,
                len(flat_values),
                c_values,
            )
            if res != 0:
                msg = f"Error in add_data_simple - error code {res} for patient_id {patient_id}, signal_name {signal_name}, data: {elem}"
                print(msg)
                messages.append(msg)
        # No return value, errors are handled by the library
        return messages

    @__test_not_disposed
    def __add_data_old_api(self, json_data: str) -> list[str]:
        """This function recieves data json object and loads the data into the algomarker patient data repository.
        Errors are collected in a string - each error in separate line. When there are no errors, the output is None.

        Notes
        -----
        The input data json request is documented in different document and the potential errors
        """

        """
        """
        js_req = json.loads(json_data)  # Check if the json is valid

        pid = int(js_req["patient_id"])
        sigs_data = js_req["signals"]
        all_data = []
        messages = []
        for sig_eme in sigs_data:
            sig_name = sig_eme["code"]
            data = sig_eme["data"]
            all_data = []
            for elem in data:
                if "timestamp" not in elem or "value" not in elem:
                    raise ValueError(
                        f"Error in data json - each signal must have 'timestamp' and 'value' fields, but got {elem}"
                    )
                timestamps = list(map(lambda x: int(x), elem["timestamp"]))
                # AddDataStr for categorical signals is not supported right now. In current algomarkersm, there are not categorical signals
                values = list(map(lambda x: float(x), elem["value"]))
                sig_data = SingleDataElement(timestamps, values)
                all_data.append(sig_data)
            res = self.add_data_simple(pid, sig_name, all_data)
            messages.extend(res)
        return messages

    @__test_not_disposed
    def add_data(self, json_data: str) -> str | None:
        """This function recieves data json object and loads the data into the algomarker patient data repository.
        Errors are collected in a string - each error in separate line. When there are no errors, the output is None.

        Notes
        -----
        The input data json request is documented in different document and the potential errors
        """
        assert self.__lib is not None
        if self.api_version == 1:
            res = self.__add_data_old_api(json_data)
            if len(res)> 0:
                return "\n".join(res)
            else:
                return None
        # For new API
        js_data = ctypes.create_string_buffer(json_data.encode("ascii"))
        res_messages = ctypes.c_char_p()
        res = self.__lib.AM_API_AddDataByType(
            self.__obj, js_data, ctypes.byref(res_messages)
        )
        if res != 0:
            print(f"AddData Failed {res}, messages ")
            res_messages_str = res_messages.value
            self.__dispose_string_mem(res_messages)
            res_messages_str_val = ""
            if res_messages_str is not None:
                res_messages_str_val = res_messages_str.decode("ascii")
            print(res_messages_str_val)
            return res_messages_str_val
        return None

    @__test_not_disposed
    def __calculate_old_api(self, request_json: str) -> dict[str, Any]:
        """Recieved json request for calculation and returns json string responde object with the result

        Notes
        -----
        The input json request and json response results are documented in a different document
        """
        assert self.__lib is not None
        # 1. Create Request Object:
        js_req = json.loads(request_json)  # Check if the json is valid
        assert (
            js_req["type"] == "request"
        )  # "Error in request json - type must be 'request'"
        request_type = ctypes.byref(ctypes.c_char_p(b"Raw"))  # Default request type
        requests = js_req["requests"]
        load_data = js_req.get("load", 0)
        pids = []
        times = []
        load_err_msgs= []
        for req in requests:
            if "patient_id" not in req or "time" not in req:
                raise ValueError(
                    "Error in request json - each request must have patient_id and time"
                )
            patient_id = int(req["patient_id"])
            time = int(req["time"])
            pids.append(patient_id)
            times.append(time)
            if load_data:
                if "data" not in req or "signals" not in req["data"]:
                    raise ValueError(
                        "Error in request json - when load is true, each request must have 'data' with 'signals'"
                    )
                load_res = self.add_data(
                    json.dumps(
                        {"signals": req["data"]["signals"], "patient_id": patient_id}
                    )
                )
                if load_res is not None:
                    load_err_msgs.extend(load_res.split("\n"))

        full_response = {
            "type": "response",
            "responses": [],
            "request_id": js_req["request_id"],
        }
        if len(load_err_msgs) > 0:
            full_response["errors"] = load_err_msgs
            return full_response
        # Convert to ctypes arrays
        c_pids = (ctypes.c_int32 * len(pids))(*pids)
        c_times = (ctypes.c_long * len(times))(*times)
        req_object = ctypes.c_void_p()
        self.__lib.AM_API_CreateRequest(
            ctypes.create_string_buffer(js_req["request_id"].encode("ascii")),
            request_type,
            1,
            c_pids,
            c_times,
            len(pids),
            ctypes.byref(req_object),
        )

        # 2. Create response object
        response_object = ctypes.c_void_p()
        self.__lib.AM_API_CreateResponses(ctypes.byref(response_object))

        # 3. Call the Calculate function
        # res_resp = ctypes.c_char_p()
        res = self.__lib.AM_API_Calculate(self.__obj, req_object, response_object)
        if res != 0:
            print(f"Error in Calculate - error code {res}")

        # 4. Check the result
        n_resp = self.__lib.AM_API_GetResponsesNum(response_object)
        print(f"Has {n_resp} responses")
        # AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
        n_msgs = ctypes.c_int32()
        msg_codes = ctypes.POINTER(ctypes.c_int32)()
        msgs_errs = ctypes.POINTER(ctypes.c_char_p)()
        res = self.__lib.AM_API_GetSharedMessages(
            response_object,
            ctypes.byref(n_msgs),  # Number of messages
            ctypes.byref(msg_codes),  # Message codes
            ctypes.byref(msgs_errs),  # Messages errors
        )
        if res != 0:
            print(f"Error in AM_API_GetSharedMessages - error code {res}")
            full_response["errors"] = [
                f"Error in AM_API_GetSharedMessages - error code {res}"
            ]
        n_msgs = n_msgs.value
        print(f"Response has {n_msgs} shared messages")
        for i in range(n_msgs):
            msg_code = msg_codes[i]
            msg_err = msgs_errs[i].decode("ascii") if msgs_errs else "None"
            if "errors" not in full_response:
                full_response["errors"] = []
            full_response["errors"].append(f"({msg_code}){msg_err}")
            print(f"Message {i}: Code: {msg_code}, Error: {msg_err}")

        for i in range(n_resp):
            # AM_API_GetResponseAtIndex(response_object, i, &response);
            # We would normally retrieve the response data here, but the old API does not provide a way to do this.
            # Here we would normally retrieve the response data, but the old API does not provide a way to do this.
            # We would need to implement the necessary functions in the C library to retrieve the response data.
            # For example:
            curr_resp_obj = ctypes.c_void_p()
            curr_num = ctypes.c_int32()
            # AM_API_GetResponseAtIndex(response_object, i, &response);
            res = self.__lib.AM_API_GetResponseAtIndex(
                response_object, i, ctypes.byref(curr_resp_obj)
            )
            if res != 0:
                print(f"Error in fetch response {i} - error code {res}")
            # AM_API_GetResponseScoresNum(response, &n_scores);
            res = self.__lib.AM_API_GetResponseScoresNum(
                curr_resp_obj, ctypes.byref(curr_num)
            )
            if res != 0:
                print(f"Error in AM_API_GetResponseScoresNum {i} - error code {res}")
            curr_num = curr_num.value
            print(f"Has {curr_num} scores in response {i}")
            # AM_API_GetResponsePoint(response, &pid, &ts);
            pid = ctypes.c_int32()
            ts = ctypes.c_long()
            res = self.__lib.AM_API_GetResponsePoint(
                curr_resp_obj, ctypes.byref(pid), ctypes.byref(ts)
            )
            if res != 0:
                print(f"Error in AM_API_GetResponsePoint {i} - error code {res}")
            pid = pid.value
            ts = ts.value
            print(f"Response {i} - Patient ID: {pid}, Timestamp: {ts}")
            n_msgs = ctypes.c_int32()
            msg_codes = ctypes.POINTER(ctypes.c_int32)()
            msgs_errs = ctypes.POINTER(ctypes.c_char_p)()
            # AM_API_GetResponseMessages(response, &n_msgs, &msg_codes, &msgs_errs);
            res = self.__lib.AM_API_GetResponseMessages(
                curr_resp_obj,
                ctypes.byref(n_msgs),  # Number of messages
                ctypes.byref(msg_codes),  # Message codes
                ctypes.byref(msgs_errs),  # Messages errors
            )
            if res != 0:
                print(f"Error in AM_API_GetResponseMessages {i} - error code {res}")
            n_msgs = n_msgs.value
            js_resp = {
                "patient_id": pid,
                "time": ts,
                "prediction": -9999,
                "messages": [],
            }

            print(f"Response {i} has {n_msgs} messages")
            for j in range(n_msgs):
                msg_code = msg_codes[j]
                msg_err = msgs_errs[j].decode("ascii") if msgs_errs else "None"
                js_resp["messages"].append(f"({msg_code}){msg_err}")
                print(f"Message {j}: Code: {msg_code}, Error: {msg_err}")
            # AM_API_GetScoreMessages
            for j in range(curr_num):
                res = self.__lib.AM_API_GetScoreMessages(
                    curr_resp_obj,
                    j,  # Assuming we want the first score messages
                    ctypes.byref(ctypes.c_int32(n_msgs)),  # Number of messages
                    ctypes.byref(msg_codes),  # Message codes
                    ctypes.byref(msgs_errs),  # Messages errors
                )
                if res != 0:
                    print(
                        f"Error in AM_API_GetScoreMessages {i} {j} - error code {res}"
                    )
            # resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
            scr_value: ctypes.c_float = ctypes.c_float()
            scr_type: ctypes.c_char_p = ctypes.c_char_p()
            for j in range(curr_num):
                res = self.__lib.AM_API_GetResponseScoreByIndex(
                    curr_resp_obj, j, ctypes.byref(scr_value), ctypes.byref(scr_type)
                )
                if res != 0:
                    print(
                        f"Error in AM_API_GetResponseScoreByIndex {i} - error code {res}"
                    )
                scr_value_v = scr_value.value
                scr_type_v = None
                if scr_type.value is not None:
                    scr_type_v = scr_type.value.decode("ascii") if scr_type else "None"
                print(
                    f"Response {i} Score {j}: Value: {scr_value_v}, Type: {scr_type_v}"
                )
                # Take the right index from exports - currently only 'pred_0' is supported for sigle pred score

                js_resp["prediction"] = scr_value_v
            full_response["responses"].append(js_resp)

        # 5. Dispose request and response objects
        self.__lib.AM_API_DisposeRequest(req_object)
        self.__lib.AM_API_DisposeResponses(response_object)
        return full_response

    @__test_not_disposed
    def calculate(self, request_json: str) -> dict[str, Any]:
        """Recieved json request for calculation and returns json string responde object with the result

        Notes
        -----
        The input json request and json response results are documented in a different document
        """
        assert self.__lib is not None
        if self.api_version == 1:
            return self.__calculate_old_api(request_json)
        js_req = ctypes.create_string_buffer(request_json.encode("ascii"))
        res_resp = ctypes.c_char_p()
        res = self.__lib.AM_API_CalculateByType(
            self.__obj, 3001, js_req, ctypes.byref(res_resp)
        )
        if res != 0:
            print(f"Calculate Failed {res}")
        try:
            res_resp_str = res_resp.value
            self.__dispose_string_mem(res_resp)
            if res_resp_str is None:
                raise NameError("Error in Calculate - response is None")
            res_resp_str = json.loads(res_resp_str)
            return res_resp_str
        except:
            print("Error in converting respond json in calculate")
            traceback.print_exc()
            raise


# Old API testing
# bdate=(ctypes.c_long * 1)(*[1988])
# bdate_right=(ctypes.c_float * 1)(*[19880327])
# am.lib.AM_API_AddData(am.obj,1,ctypes.create_string_buffer(b"BDATE"),1, bdate,0 ,ctypes.POINTER(ctypes.c_float)())
# am.lib.AM_API_AddData(am.obj,1,ctypes.create_string_buffer(b"BDATE"),0, ctypes.POINTER(ctypes.c_long)(),1 ,bdate_right)

if __name__ == "__main__":
    print(
        "This is a module for AlgoMarker Python API. Use it as a module, not as a script."
    )
    print("Example usage:")
    AlgoMarker_path = os.path.join(
        os.environ["HOME"],
        "Documents/MES/AlgoMarkers/AM_LGI/AlgoMarker/ColonFlag_3.1.0.0/ColonFlag-3.1.amconfig",
        # "Documents/MES/AlgoMarkers/docker_images/LGI-Flag-ButWhy-3.1.2-Scorer/data/app/LGI-Flag-ButWhy-3.1.2-Scorer/LGI-ColonFlag-3.1.amconfig"
    )
    libpath = None
    libpath = os.path.join(
        os.environ["HOME"],
        "Documents/MES/AlgoMarkers/AM_LGI/AlgoMarker/ColonFlag_3.1.0.0/libdyn_AlgoMarker.25102018_1.so",
        # "Documents/MES/AlgoMarkers/docker_images/LGI-Flag-ButWhy-3.1.2-Scorer/data/app/LGI-Flag-ButWhy-3.1.2-Scorer/lib/libdyn_AlgoMarker.so"
    )
    request_json = AlgoMarker.create_request_json(1, 20240101)
    with AlgoMarker(AlgoMarker_path, libpath) as am:
        print(am.discovery())
        am.clear_data()
        am.add_data_simple(1, "BYEAR", [SingleDataElement([], [1978])])
        am.add_data_simple(1, "GENDER", [SingleDataElement([], [1])])
        am.add_data_simple(
            1,
            "Hemoglobin",
            [
                SingleDataElement([20220101], [14.5]),
                SingleDataElement([20230101], [14.5]),
                SingleDataElement([20240101], [14.5]),
            ],
        )
        am.add_data_simple(
            1,
            "Hematocrit",
            [
                SingleDataElement([20220101], [33]),
                SingleDataElement([20230101], [33]),
                SingleDataElement([20240101], [33]),
            ],
        )
        am.add_data_simple(
            1,
            "MCH",
            [
                SingleDataElement([20220101], [33]),
                SingleDataElement([20230101], [33]),
                SingleDataElement([20240101], [33]),
            ],
        )
        am.add_data_simple(
            1,
            "RBC",
            [
                SingleDataElement([20220101], [4.5]),
                SingleDataElement([20230101], [4.5]),
                SingleDataElement([20240101], [4.5]),
            ],
        )
        am.add_data_simple(
            1,
            "MCV",
            [
                SingleDataElement([20220101], [90]),
                SingleDataElement([20230101], [90]),
                SingleDataElement([20240101], [90]),
            ],
        )
        resp = am.calculate(request_json)
        print("Response:")
        print(resp)
        print("Done with AlgoMarker example")
