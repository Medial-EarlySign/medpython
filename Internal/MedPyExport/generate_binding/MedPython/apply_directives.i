%rename(Time) MPTime;
%rename(Features) MPFeatures;
%rename(Sample) MPSample;
%rename(SampleVectorAdaptor) MPSampleVectorAdaptor;
%rename(IdSamples) MPIdSamples;
%rename(IdSamplesVectorAdaptor) MPIdSamplesVectorAdaptor;
%rename(Samples) MPSamples;
%apply (int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** ids, int* num_ids)}
%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float** out_predbuf, int* out_predbuf_len), (float** preds_buf, int* preds_buf_len), (float** y_buf, int* y_buf_len), (float** categs_buf, int* categs_buf_len)}
%apply (float* IN_ARRAY1, int DIM1) {(float* in_predbuf, int in_predbuf_len)}
%extend MPIdSamples{
  %pythoncode %{
    __swig_getmethods__["id"] = MEDPY_GET_id
    __swig_setmethods__["id"] = MEDPY_SET_id
    __swig_getmethods__["split"] = MEDPY_GET_split
    __swig_setmethods__["split"] = MEDPY_SET_split
    __swig_getmethods__["samples"] = MEDPY_GET_samples
    if _newclass:
      id = property(MEDPY_GET_id,MEDPY_SET_id)
      samples = property(MEDPY_GET_samples,None)
      split = property(MEDPY_GET_split,MEDPY_SET_split)
  %}
};
%extend MPSample{
  %pythoncode %{
    __swig_getmethods__["id"] = MEDPY_GET_id
    __swig_setmethods__["id"] = MEDPY_SET_id
    __swig_getmethods__["split"] = MEDPY_GET_split
    __swig_setmethods__["split"] = MEDPY_SET_split
    __swig_getmethods__["time"] = MEDPY_GET_time
    __swig_setmethods__["time"] = MEDPY_SET_time
    __swig_getmethods__["outcome"] = MEDPY_GET_outcome
    __swig_setmethods__["outcome"] = MEDPY_SET_outcome
    __swig_getmethods__["outcomeTime"] = MEDPY_GET_outcomeTime
    __swig_setmethods__["outcomeTime"] = MEDPY_SET_outcomeTime
    __swig_getmethods__["prediction"] = MEDPY_GET_prediction
    __swig_setmethods__["prediction"] = MEDPY_SET_prediction
    if _newclass:
      outcomeTime = property(MEDPY_GET_outcomeTime,MEDPY_SET_outcomeTime)
      prediction = property(MEDPY_GET_prediction,MEDPY_SET_prediction)
      split = property(MEDPY_GET_split,MEDPY_SET_split)
      time = property(MEDPY_GET_time,MEDPY_SET_time)
      outcome = property(MEDPY_GET_outcome,MEDPY_SET_outcome)
      id = property(MEDPY_GET_id,MEDPY_SET_id)
  %}
};
%extend MPSamples{
  %pythoncode %{
    __swig_getmethods__["time_unit"] = MEDPY_GET_time_unit
    __swig_setmethods__["time_unit"] = MEDPY_SET_time_unit
    __swig_getmethods__["idSamples"] = MEDPY_GET_idSamples
    if _newclass:
      time_unit = property(MEDPY_GET_time_unit,MEDPY_SET_time_unit)
      idSamples = property(MEDPY_GET_idSamples,None)
  %}
};
%rename(PidRepository) MPPidRepository;
%rename(Dictionary) MPDictionary;
%apply (char** ARGOUTVIEWM_ARRAY1, int* DIM1) {(char** lut_array, int* lut_size)}
%apply (int* IN_ARRAY1, int DIM1) {(int* members_array, int members_size)}
%rename(FeatureAttr) MPFeatureAttr;
%rename(StringFeatureAttrMapAdaptor) MPStringFeatureAttrMapAdaptor;
%apply (int* IN_ARRAY1, int DIM1) {(int* int_in_buf, int int_in_buf_len)}
%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float** float_out_buf, int* float_out_buf_len)}
%extend MPFeatureAttr{
  %pythoncode %{
    __swig_setmethods__["normalized"] = MEDPY_SET_normalized
    __swig_getmethods__["normalized"] = MEDPY_GET_normalized
    __swig_setmethods__["imputed"] = MEDPY_SET_imputed
    __swig_getmethods__["imputed"] = MEDPY_GET_imputed
    if _newclass:
      normalized = property(MEDPY_GET_normalized,MEDPY_SET_normalized)
      imputed = property(MEDPY_GET_imputed,MEDPY_SET_imputed)
  %}
};
%extend MPFeatures{
  %pythoncode %{
    __swig_getmethods__["data"] = MEDPY_GET_data
    __swig_getmethods__["weights"] = MEDPY_GET_weights
    __swig_getmethods__["samples"] = MEDPY_GET_samples
    __swig_getmethods__["pid_pos_len"] = MEDPY_GET_pid_pos_len
    __swig_getmethods__["attributes"] = MEDPY_GET_attributes
    __swig_getmethods__["tags"] = MEDPY_GET_tags
    __swig_getmethods__["time_unit"] = MEDPY_GET_time_unit
    __swig_setmethods__["time_unit"] = MEDPY_SET_time_unit
    __swig_getmethods__["global_serial_id_cnt"] = MEDPY_GET_global_serial_id_cnt
    __swig_setmethods__["global_serial_id_cnt"] = MEDPY_SET_global_serial_id_cnt
    if _newclass:
      global_serial_id_cnt = property(MEDPY_GET_global_serial_id_cnt,MEDPY_SET_global_serial_id_cnt)
      tags = property(MEDPY_GET_tags,None)
      pid_pos_len = property(MEDPY_GET_pid_pos_len,None)
      weights = property(MEDPY_GET_weights,None)
      time_unit = property(MEDPY_GET_time_unit,MEDPY_SET_time_unit)
      samples = property(MEDPY_GET_samples,None)
      attributes = property(MEDPY_GET_attributes,None)
      data = property(MEDPY_GET_data,None)
  %}
};
%rename(SigVectorAdaptor) MPSigVectorAdaptor;
%rename(Sig) MPSig;
%rename(SigExporter) MPSigExporter;
%apply (int* IN_ARRAY1, int DIM1) {(int* pids_to_take, int num_pids_to_take), (int* pids_to_take, int num_pids_to_take)}
%extend MPPidRepository{
  %pythoncode %{
    __swig_getmethods__["pids"] = MEDPY_GET_pids
    if _newclass:
      pids = property(MEDPY_GET_pids,None)
  %}
};
%extend MPSigVectorAdaptor{
  %pythoncode %{
    __swig_getmethods__["type"] = MEDPY_GET_type
    __swig_getmethods__["n_time_channels"] = MEDPY_GET_n_time_channels
    __swig_getmethods__["n_val_channels"] = MEDPY_GET_n_val_channels
    __swig_getmethods__["time_unit"] = MEDPY_GET_time_unit
    __swig_getmethods__["size"] = MEDPY_GET_size
    if _newclass:
      size = property(MEDPY_GET_size,None)
      type = property(MEDPY_GET_type,None)
      time_unit = property(MEDPY_GET_time_unit,None)
      n_time_channels = property(MEDPY_GET_n_time_channels,None)
      n_val_channels = property(MEDPY_GET_n_val_channels,None)
  %}
};
%rename(IntIntMapAdaptor) MPIntIntMapAdaptor;
%rename(StringStringMapAdaptor) MPStringStringMapAdaptor;
%rename(StringVecFloatMapAdaptor) MPStringVecFloatMapAdaptor;
%rename(IntPairIntIntMapAdaptor) MPIntPairIntIntMapAdaptor;
%rename(StringUOSetStringMapAdaptor) MPStringUOSetStringMapAdaptor;
%rename(IntStringMapAdaptor) MPIntStringMapAdaptor;
%rename(IntVecIntMapAdaptor) MPIntVecIntMapAdaptor;
%apply (int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** intkeys_out_buf, int* intkeys_out_buf_len), (int** int_out_buf, int* int_out_buf_len), (int** intkeys_out_buf, int* intkeys_out_buf_len), (int** int_out_buf, int* int_out_buf_len)}
%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float** float_out_buf, int* float_out_buf_len)}
%apply (float* IN_ARRAY1, int DIM1) {(float* float_in_buf, int float_in_buf_len)}
%apply (int* IN_ARRAY1, int DIM1) {(int* int_in_buf, int int_in_buf_len), (int* int_in_buf, int int_in_buf_len)}
%rename(ModelStage) MPModelStage;
%rename(Model) MPModel;
%extend MPModel{
  %pythoncode %{
    __swig_getmethods__["features"] = MEDPY_GET_features
    __swig_getmethods__["verbosity"] = MEDPY_GET_verbosity
    __swig_setmethods__["verbosity"] = MEDPY_SET_verbosity
    if _newclass:
      verbosity = property(MEDPY_GET_verbosity,MEDPY_SET_verbosity)
      features = property(MEDPY_GET_features,None)
  %}
};
%rename(Split) MPSplit;
%extend MPSplit{
  %pythoncode %{
    __swig_getmethods__["nsplits"] = MEDPY_GET_nsplits
    __swig_getmethods__["pid2split"] = MEDPY_GET_pid2split
    if _newclass:
      pid2split = property(MEDPY_GET_pid2split,None)
      nsplits = property(MEDPY_GET_nsplits,None)
  %}
};
%rename(Mat) MPMat;
%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {(float** rowv, int* rowv_n), (float** colv, int* colv_n), (float** buf_avg, int* buf_avg_n), (float** buf_std, int* buf_std_n)}
%apply (int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** out_row_ids_buf, int* out_row_ids_buf_len), (int** avg_buf, int* avg_buf_len), (int** std_buf, int* std_buf_len)}
%apply (float* IN_ARRAY1, int DIM1) {(float* m_add, int nrows_to_add), (float* m_add, int ncols_to_add), (float* wgts, int wgts_n), (float* external_avg, int external_avg_n), (float* external_std, int external_std_n)}
%apply (int* IN_ARRAY1, int DIM1) {(int* row_ids_buf, int row_ids_buf_len)}
%extend MPMat{
  %pythoncode %{
    __swig_getmethods__["Normalize_Cols"] = MEDPY_GET_Normalize_Cols
    __swig_getmethods__["Normalize_Rows"] = MEDPY_GET_Normalize_Rows
    __swig_getmethods__["nrows"] = MEDPY_GET_nrows
    __swig_getmethods__["ncols"] = MEDPY_GET_ncols
    __swig_getmethods__["row_ids"] = MEDPY_GET_row_ids
    __swig_setmethods__["row_ids"] = MEDPY_SET_row_ids
    __swig_getmethods__["signals"] = MEDPY_GET_signals
    __swig_getmethods__["avg"] = MEDPY_GET_avg
    __swig_getmethods__["std"] = MEDPY_GET_std
    __swig_getmethods__["normalized_flag"] = MEDPY_GET_normalized_flag
    __swig_setmethods__["normalized_flag"] = MEDPY_SET_normalized_flag
    __swig_getmethods__["transposed_flag"] = MEDPY_GET_transposed_flag
    __swig_setmethods__["transposed_flag"] = MEDPY_SET_transposed_flag
    __swig_getmethods__["missing_value"] = MEDPY_GET_missing_value
    __swig_setmethods__["missing_value"] = MEDPY_SET_missing_value
    if _newclass:
      row_ids = property(MEDPY_GET_row_ids,MEDPY_SET_row_ids)
      signals = property(MEDPY_GET_signals,None)
      Normalize_Cols = property(MEDPY_GET_Normalize_Cols,None)
      Normalize_Rows = property(MEDPY_GET_Normalize_Rows,None)
      ncols = property(MEDPY_GET_ncols,None)
      normalized_flag = property(MEDPY_GET_normalized_flag,MEDPY_SET_normalized_flag)
      missing_value = property(MEDPY_GET_missing_value,MEDPY_SET_missing_value)
      nrows = property(MEDPY_GET_nrows,None)
      std = property(MEDPY_GET_std,None)
      transposed_flag = property(MEDPY_GET_transposed_flag,MEDPY_SET_transposed_flag)
      avg = property(MEDPY_GET_avg,None)
  %}
};
%apply (int* IN_ARRAY1, int DIM1) {(int* pids_to_take, int num_pids_to_take)}
%apply (void** ARGOUTMVAR_ARRAY1, int* DIM1, int* NPYDTC1) {(void** outarr1, int* outarr1_sz, int* outarr1_npytype), (void** outarr1, int* outarr1_sz, int* outarr1_npytype), (void** outarr1, int* outarr1_sz, int* outarr1_npytype)}
%rename(PandasAdaptor) MPPandasAdaptor;
%apply (void** ARGOUTMVAR_ARRAY1, int* DIM1, int* NPYDTC1) {(void** outarr1, int* outarr1_sz, int* outarr1_npytype), (void** outarr1, int* outarr1_sz, int* outarr1_npytype), (void** outarr1, int* outarr1_sz, int* outarr1_npytype)}
