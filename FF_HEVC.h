#ifndef FF_HEVC_H
#define FF_HEVC_H
enum MS_HEVCNAL_UNITTYPE
{
	MS_HEVC_NAL_TRAIL_N = 0,
	MS_HEVC_NAL_TRAIL_R = 1,
	MS_HEVC_NAL_TSA_N = 2,
	MS_HEVC_NAL_TSA_R = 3,
	MS_HEVC_NAL_STSA_N = 4,
	MS_HEVC_NAL_STSA_R = 5,
	MS_HEVC_NAL_RADL_N = 6,
	MS_HEVC_NAL_RADL_R = 7,
	MS_HEVC_NAL_RASL_N = 8,
	MS_HEVC_NAL_RASL_R = 9,
	MS_HEVC_NAL_VCL_N10 = 10,
	MS_HEVC_NAL_VCL_R11 = 11,
	MS_HEVC_NAL_VCL_N12 = 12,
	MS_HEVC_NAL_VCL_R13 = 13,
	MS_HEVC_NAL_VCL_N14 = 14,
	MS_HEVC_NAL_VCL_R15 = 15,
	MS_HEVC_NAL_BLA_W_LP = 16,
	MS_HEVC_NAL_BLA_W_RADL = 17,
	MS_HEVC_NAL_BLA_N_LP = 18,
	MS_HEVC_NAL_IDR_W_RADL = 19,
	MS_HEVC_NAL_IDR_N_LP = 20,
	MS_HEVC_NAL_CRA_NUT = 21,
	MS_HEVC_NAL_IRAP_VCL22 = 22,
	MS_HEVC_NAL_IRAP_VCL23 = 23,
	MS_HEVC_NAL_RSV_VCL24 = 24,
	MS_HEVC_NAL_RSV_VCL25 = 25,
	MS_HEVC_NAL_RSV_VCL26 = 26,
	MS_HEVC_NAL_RSV_VCL27 = 27,
	MS_HEVC_NAL_RSV_VCL28 = 28,
	MS_HEVC_NAL_RSV_VCL29 = 29,
	MS_HEVC_NAL_RSV_VCL30 = 30,
	MS_HEVC_NAL_RSV_VCL31 = 31,
	MS_HEVC_NAL_VPS = 32,
	MS_HEVC_NAL_SPS = 33,
	MS_HEVC_NAL_PPS = 34,
	MS_HEVC_NAL_AUD = 35,
	MS_HEVC_NAL_EOS_NUT = 36,
	MS_HEVC_NAL_EOB_NUT = 37,
	MS_HEVC_NAL_FD_NUT = 38,
	MS_HEVC_NAL_SEI_PREFIX = 39,
	MS_HEVC_NAL_SEI_SUFFIX = 40,
};
enum
{
	MS_HEVC_MAX_LAYERS = 63,
	MS_HEVC_MAX_SUB_LAYERS = 7,
	MS_HEVC_MAX_LAYER_SETS = 1024,
	MS_HEVC_MAX_VPS_COUNT = 16,
	MS_HEVC_MAX_SPS_COUNT = 16,
	MS_HEVC_MAX_PPS_COUNT = 64,
	MS_HEVC_MAX_DPB_SIZE = 16,
	MS_HEVC_MAX_REFS = MS_HEVC_MAX_DPB_SIZE,
	MS_HEVC_MAX_SHORT_TERM_REF_PIC_SETS = 64,
	MS_HEVC_MAX_LONG_TERM_REF_PICS = 32,
	MS_HEVC_MIN_LOG2_CTB_SIZE = 4,
	MS_HEVC_MAX_LOG2_CTB_SIZE = 6,
	MS_HEVC_MAX_CPB_CNT = 32,
	MS_HEVC_MAX_LUMA_PS = 35651584,
	MS_HEVC_MAX_WIDTH = 16888,
	MS_HEVC_MAX_HEIGHT = 16888,
	MS_HEVC_MAX_TILE_ROWS = 22,
	MS_HEVC_MAX_TILE_COLUMNS = 20,
	MS_HEVC_MAX_ENTRY_POINT_OFFSETS = MS_HEVC_MAX_TILE_COLUMNS * 135,
};
typedef struct HVCCNALUnitArray
{
	UINT8	array_completeness;
	UINT8	NAL_unit_type;
	UINT16	numNalus;
	UINT16* nalUnitLength;
	UINT8** nalUnit;
} HVCCNALUnitArray;

typedef struct HEVCDecoderConfigurationRecord
{
	UINT8  configurationVersion;
	UINT8  general_profile_space;
	UINT8  general_tier_flag;
	UINT8  general_profile_idc;
	UINT32 general_profile_compatibility_flags;
	UINT64 general_constraint_indicator_flags;
	UINT8  general_level_idc;
	UINT16 min_spatial_segmentation_idc;
	UINT8  parallelismType;
	UINT8  chromaFormat;
	UINT8  bitDepthLumaMinus8;
	UINT8  bitDepthChromaMinus8;
	UINT16 avgFrameRate;
	UINT8  constantFrameRate;
	UINT8  numTemporalLayers;
	UINT8  temporalIdNested;
	UINT8  lengthSizeMinusOne;
	UINT8  numOfArrays;
	HVCCNALUnitArray* array;
} HEVCDecoderConfigurationRecord;

typedef struct HVCCProfileTierLevel
{
	UINT8  profile_space;
	UINT8  tier_flag;
	UINT8  profile_idc;
	UINT32 profile_compatibility_flags;
	UINT64 constraint_indicator_flags;
	UINT8  level_idc;
} HVCCProfileTierLevel;

enum
{
	HEVC_MAX_LAYERS = 63,
	HEVC_MAX_SUB_LAYERS = 7,
	HEVC_MAX_LAYER_SETS = 1024,
	HEVC_MAX_VPS_COUNT = 16,
	HEVC_MAX_SPS_COUNT = 16,
	HEVC_MAX_PPS_COUNT = 64,
	HEVC_MAX_DPB_SIZE = 16,
	HEVC_MAX_REFS = HEVC_MAX_DPB_SIZE,
	HEVC_MAX_SHORT_TERM_REF_PIC_SETS = 64,
	HEVC_MAX_LONG_TERM_REF_PICS = 32,
	HEVC_MIN_LOG2_CTB_SIZE = 4,
	HEVC_MAX_LOG2_CTB_SIZE = 6,
	HEVC_MAX_CPB_CNT = 32,
	HEVC_MAX_LUMA_PS = 35651584,
	HEVC_MAX_WIDTH = 16888,
	HEVC_MAX_HEIGHT = 16888,
	HEVC_MAX_TILE_ROWS = 22,
	HEVC_MAX_TILE_COLUMNS = 20,
	HEVC_MAX_SLICE_SEGMENTS = 600,
	HEVC_MAX_ENTRY_POINT_OFFSETS = HEVC_MAX_TILE_COLUMNS * 135,
};

static const uint8_t hevc_sub_width_c[] =
{
	1, 2, 2, 1
};

static const uint8_t hevc_sub_height_c[] =
{
	1, 2, 1, 1
};

typedef struct ShortTermRPS
{
	unsigned int num_negative_pics;
	int num_delta_pocs;
	int rps_idx_num_delta_pocs;
	int32_t delta_poc[32];
	uint8_t used[32];
} ShortTermRPS;

typedef struct HEVCWindow
{
	unsigned int left_offset;
	unsigned int right_offset;
	unsigned int top_offset;
	unsigned int bottom_offset;
} HEVCWindow;

typedef struct VUI
{
	AVRational sar;
	int overscan_info_present_flag;
	int overscan_appropriate_flag;
	int video_signal_type_present_flag;
	int video_format;
	int video_full_range_flag;
	int colour_description_present_flag;
	uint8_t colour_primaries;
	uint8_t transfer_characteristic;
	uint8_t matrix_coeffs;
	int chroma_loc_info_present_flag;
	int chroma_sample_loc_type_top_field;
	int chroma_sample_loc_type_bottom_field;
	int neutra_chroma_indication_flag;
	int field_seq_flag;
	int frame_field_info_present_flag;
	int default_display_window_flag;
	HEVCWindow def_disp_win;
	int vui_timing_info_present_flag;
	uint32_t vui_num_units_in_tick;
	uint32_t vui_time_scale;
	int vui_poc_proportional_to_timing_flag;
	int vui_num_ticks_poc_diff_one_minus1;
	int vui_hrd_parameters_present_flag;
	int bitstream_restriction_flag;
	int tiles_fixed_structure_flag;
	int motion_vectors_over_pic_boundaries_flag;
	int restricted_ref_pic_lists_flag;
	int min_spatial_segmentation_idc;
	int max_bytes_per_pic_denom;
	int max_bits_per_min_cu_denom;
	int log2_max_mv_length_horizontal;
	int log2_max_mv_length_vertical;
} VUI;

typedef struct PTLCommon
{
	uint8_t profile_space;
	uint8_t tier_flag;
	uint8_t profile_idc;
	uint8_t profile_compatibility_flag[32];
	uint8_t progressive_source_flag;
	uint8_t interlaced_source_flag;
	uint8_t non_packed_constraint_flag;
	uint8_t frame_only_constraint_flag;
	uint8_t max_12bit_constraint_flag;
	uint8_t max_10bit_constraint_flag;
	uint8_t max_8bit_constraint_flag;
	uint8_t max_422chroma_constraint_flag;
	uint8_t max_420chroma_constraint_flag;
	uint8_t max_monochrome_constraint_flag;
	uint8_t intra_constraint_flag;
	uint8_t one_picture_only_constraint_flag;
	uint8_t lower_bit_rate_constraint_flag;
	uint8_t max_14bit_constraint_flag;
	uint8_t inbld_flag;
	uint8_t level_idc;
} PTLCommon;

typedef struct PTL
{
	PTLCommon general_ptl;
	PTLCommon sub_layer_ptl[HEVC_MAX_SUB_LAYERS];
	uint8_t sub_layer_profile_present_flag[HEVC_MAX_SUB_LAYERS];
	uint8_t sub_layer_level_present_flag[HEVC_MAX_SUB_LAYERS];
} PTL;

typedef struct tagHEVCVPS
{
	uint8_t vps_temporal_id_nesting_flag;
	int vps_max_layers;
	int vps_max_sub_layers; ///< vps_max_temporal_layers_minus1 + 1
	PTL ptl;
	int vps_sub_layer_ordering_info_present_flag;
	unsigned int vps_max_dec_pic_buffering[HEVC_MAX_SUB_LAYERS];
	unsigned int vps_num_reorder_pics[HEVC_MAX_SUB_LAYERS];
	unsigned int vps_max_latency_increase[HEVC_MAX_SUB_LAYERS];
	int vps_max_layer_id;
	int vps_num_layer_sets; ///< vps_num_layer_sets_minus1 + 1
	uint8_t vps_timing_info_present_flag;
	uint32_t vps_num_units_in_tick;
	uint32_t vps_time_scale;
	uint8_t vps_poc_proportional_to_timing_flag;
	int vps_num_ticks_poc_diff_one; ///< vps_num_ticks_poc_diff_one_minus1 + 1
	int vps_num_hrd_parameters;
	uint8_t data[4096];
	int data_size;
} HEVCVPS;

typedef struct tagLongTermRPS
{
	int     poc[32];
	uint8_t poc_msb_present[32];
	uint8_t used[32];
	uint8_t nb_refs;
} LongTermRPS;

typedef struct tagScalingList
{
	uint8_t sl[4][6][64];
	uint8_t sl_dc[2][6];
} ScalingList;

typedef struct tagHEVCSPS
{
	unsigned vps_id;
	int chroma_format_idc;
	uint8_t separate_colour_plane_flag;
	HEVCWindow output_window;
	HEVCWindow pic_conf_win;
	int bit_depth;
	int bit_depth_chroma;
	int pixel_shift;
	enum AVPixelFormat pix_fmt;
	unsigned int log2_max_poc_lsb;
	int pcm_enabled_flag;
	int max_sub_layers;
	struct
	{
		int max_dec_pic_buffering;
		int num_reorder_pics;
		int max_latency_increase;
	} temporal_layer[HEVC_MAX_SUB_LAYERS];
	uint8_t temporal_id_nesting_flag;
	VUI vui;
	PTL ptl;
	uint8_t scaling_list_enable_flag;
	ScalingList scaling_list;
	unsigned int nb_st_rps;
	ShortTermRPS st_rps[HEVC_MAX_SHORT_TERM_REF_PIC_SETS];
	uint8_t amp_enabled_flag;
	uint8_t sao_enabled;
	uint8_t long_term_ref_pics_present_flag;
	uint16_t lt_ref_pic_poc_lsb_sps[HEVC_MAX_LONG_TERM_REF_PICS];
	uint8_t used_by_curr_pic_lt_sps_flag[HEVC_MAX_LONG_TERM_REF_PICS];
	uint8_t num_long_term_ref_pics_sps;
	struct
	{
		uint8_t bit_depth;
		uint8_t bit_depth_chroma;
		unsigned int log2_min_pcm_cb_size;
		unsigned int log2_max_pcm_cb_size;
		uint8_t loop_filter_disable_flag;
	} pcm;
	uint8_t sps_temporal_mvp_enabled_flag;
	uint8_t sps_strong_intra_smoothing_enable_flag;
	unsigned int log2_min_cb_size;
	unsigned int log2_diff_max_min_coding_block_size;
	unsigned int log2_min_tb_size;
	unsigned int log2_max_trafo_size;
	unsigned int log2_ctb_size;
	unsigned int log2_min_pu_size;
	int max_transform_hierarchy_depth_inter;
	int max_transform_hierarchy_depth_intra;
	int sps_range_extension_flag;
	int transform_skip_rotation_enabled_flag;
	int transform_skip_context_enabled_flag;
	int implicit_rdpcm_enabled_flag;
	int explicit_rdpcm_enabled_flag;
	int extended_precision_processing_flag;
	int intra_smoothing_disabled_flag;
	int high_precision_offsets_enabled_flag;
	int persistent_rice_adaptation_enabled_flag;
	int cabac_bypass_alignment_enabled_flag;
	int width;
	int height;
	int ctb_width;
	int ctb_height;
	int ctb_size;
	int min_cb_width;
	int min_cb_height;
	int min_tb_width;
	int min_tb_height;
	int min_pu_width;
	int min_pu_height;
	int tb_mask;
	int hshift[3];
	int vshift[3];
	int qp_bd_offset;
	uint8_t data[4096];
	int data_size;
} HEVCSPS;

typedef struct HEVCPPS
{
	unsigned int sps_id;
	uint8_t sign_data_hiding_flag;
	uint8_t cabac_init_present_flag;
	int num_ref_idx_l0_default_active;
	int num_ref_idx_l1_default_active;
	int pic_init_qp_minus26;
	uint8_t constrained_intra_pred_flag;
	uint8_t transform_skip_enabled_flag;
	uint8_t cu_qp_delta_enabled_flag;
	int diff_cu_qp_delta_depth;
	int cb_qp_offset;
	int cr_qp_offset;
	uint8_t pic_slice_level_chroma_qp_offsets_present_flag;
	uint8_t weighted_pred_flag;
	uint8_t weighted_bipred_flag;
	uint8_t output_flag_present_flag;
	uint8_t transquant_bypass_enable_flag;
	uint8_t dependent_slice_segments_enabled_flag;
	uint8_t tiles_enabled_flag;
	uint8_t entropy_coding_sync_enabled_flag;
	uint16_t num_tile_columns;   ///< num_tile_columns_minus1 + 1
	uint16_t num_tile_rows;      ///< num_tile_rows_minus1 + 1
	uint8_t uniform_spacing_flag;
	uint8_t loop_filter_across_tiles_enabled_flag;
	uint8_t seq_loop_filter_across_slices_enabled_flag;
	uint8_t deblocking_filter_control_present_flag;
	uint8_t deblocking_filter_override_enabled_flag;
	uint8_t disable_dbf;
	int beta_offset;    ///< beta_offset_div2 * 2
	int tc_offset;      ///< tc_offset_div2 * 2
	uint8_t scaling_list_data_present_flag;
	ScalingList scaling_list;
	uint8_t lists_modification_present_flag;
	int log2_parallel_merge_level; ///< log2_parallel_merge_level_minus2 + 2
	int num_extra_slice_header_bits;
	uint8_t slice_header_extension_present_flag;
	uint8_t log2_max_transform_skip_block_size;
	uint8_t pps_range_extensions_flag;
	uint8_t cross_component_prediction_enabled_flag;
	uint8_t chroma_qp_offset_list_enabled_flag;
	uint8_t diff_cu_chroma_qp_offset_depth;
	uint8_t chroma_qp_offset_list_len_minus1;
	int8_t  cb_qp_offset_list[6];
	int8_t  cr_qp_offset_list[6];
	uint8_t log2_sao_offset_scale_luma;
	uint8_t log2_sao_offset_scale_chroma;
	// Inferred parameters
	unsigned int* column_width;  ///< ColumnWidth
	unsigned int* row_height;    ///< RowHeight
	unsigned int* col_bd;        ///< ColBd
	unsigned int* row_bd;        ///< RowBd
	int* col_idxX;
	int* ctb_addr_rs_to_ts; ///< CtbAddrRSToTS
	int* ctb_addr_ts_to_rs; ///< CtbAddrTSToRS
	int* tile_id;           ///< TileId
	int* tile_pos_rs;       ///< TilePosRS
	int* min_tb_addr_zs;    ///< MinTbAddrZS
	int* min_tb_addr_zs_tab;///< MinTbAddrZS
	uint8_t data[4096];
	int data_size;
} HEVCPPS;

typedef struct HEVCParamSets
{
	AVBufferRef* vps_list[HEVC_MAX_VPS_COUNT];
	AVBufferRef* sps_list[HEVC_MAX_SPS_COUNT];
	AVBufferRef* pps_list[HEVC_MAX_PPS_COUNT];
	const HEVCVPS* vps;
	const HEVCSPS* sps;
	const HEVCPPS* pps;
} HEVCParamSets;

enum HEVCSliceType
{
	HEVC_SLICE_B = 0,
	HEVC_SLICE_P = 1,
	HEVC_SLICE_I = 2,
	HEVC_SLICE_NONE = 3
};

typedef struct tagHEVC_SLICE_HERDER
{
	UINT32 pps_id;
	UINT32 slice_segment_addr;
	UINT32 slice_addr;
	enum HEVCSliceType slice_type;
	INT32 pic_order_cnt_lsb;
	UINT8 first_slice_in_pic_flag;
	UINT8 dependent_slice_segment_flag;
	UINT8 pic_output_flag;
	UINT8 colour_plane_id;
	INT32 short_term_ref_pic_set_sps_flag;
	INT32 short_term_ref_pic_set_size;
	ShortTermRPS slice_rps;
	const ShortTermRPS* short_term_rps;
	INT32 long_term_ref_pic_set_size;
	LongTermRPS long_term_rps;
	UINT32 list_entry_lx[2][32];
	UINT8 rpl_modification_flag[2];
	UINT8 no_output_of_prior_pics_flag;
	UINT8 slice_temporal_mvp_enabled_flag;
	UINT32 nb_refs[2];
	UINT8 slice_sample_adaptive_offset_flag[3];
	UINT8 mvd_l1_zero_flag;
	UINT8 cabac_init_flag;
	UINT8 disable_deblocking_filter_flag;
	UINT8 slice_loop_filter_across_slices_enabled_flag;
	UINT8 collocated_list;
	UINT32 collocated_ref_idx;
	INT32 slice_qp_delta;
	INT32 slice_cb_qp_offset;
	INT32 slice_cr_qp_offset;
	UINT8 cu_chroma_qp_offset_enabled_flag;
	INT32 beta_offset;
	INT32 tc_offset;
	UINT32 max_num_merge_cand;
	unsigned* entry_point_offset;
	INT32* offset;
	INT32* size;
	INT32 num_entry_point_offsets;
	INT8 slice_qp;
	UINT8 luma_log2_weight_denom;
	INT16 chroma_log2_weight_denom;
	INT16 luma_weight_l0[16];
	INT16 chroma_weight_l0[16][2];
	INT16 chroma_weight_l1[16][2];
	INT16 luma_weight_l1[16];
	INT16 luma_offset_l0[16];
	INT16 chroma_offset_l0[16][2];
	INT16 luma_offset_l1[16];
	INT16 chroma_offset_l1[16][2];
	INT32 slice_ctb_addr_rs;
} HEVC_SLICE_HERDER;


#define MAX_SPATIAL_SEGMENTATION 4096 
#define IS_IDR(s) ((s)->nal_unit_type == MS_HEVC_NAL_IDR_W_RADL || (s)->nal_unit_type == MS_HEVC_NAL_IDR_N_LP)
#define IS_BLA(s) ((s)->nal_unit_type == MS_HEVC_NAL_BLA_W_RADL || (s)->nal_unit_type == MS_HEVC_NAL_BLA_W_LP || \
                   (s)->nal_unit_type == MS_HEVC_NAL_BLA_N_LP)
#define IS_IRAP(s) ((s)->nal_unit_type >= 16 && (s)->nal_unit_type <= 23)

static int hvcc_array_add_nal_unit(uint8_t *nal_buf, uint32_t nal_size, uint8_t nal_type, int ps_array_completeness, HEVCDecoderConfigurationRecord *hvcc)
{
	int ret;
	uint8_t index;
	uint16_t numNalus;
	HVCCNALUnitArray *array;
	for (index = 0; index < hvcc->numOfArrays; index++)
		if (hvcc->array[index].NAL_unit_type == nal_type)
			break;
	if (index >= hvcc->numOfArrays)
	{
		uint8_t i;
		ret = av_reallocp_array(&hvcc->array, index + 1, sizeof(HVCCNALUnitArray));
		if (ret < 0)
			return ret;
		for (i = hvcc->numOfArrays; i <= index; i++)
			memset(&hvcc->array[i], 0, sizeof(HVCCNALUnitArray));
		hvcc->numOfArrays = index + 1;
	}
	array = &hvcc->array[index];
	numNalus = array->numNalus;
	ret = av_reallocp_array(&array->nalUnit, numNalus + 1, sizeof(uint8_t*));
	if (ret < 0)
		return ret;
	ret = av_reallocp_array(&array->nalUnitLength, numNalus + 1, sizeof(uint16_t));
	if (ret < 0)
		return ret;
	array->nalUnit[numNalus] = nal_buf;
	array->nalUnitLength[numNalus] = nal_size;
	array->NAL_unit_type = nal_type;
	array->numNalus++;
	if (nal_type == MS_HEVC_NAL_VPS || nal_type == MS_HEVC_NAL_SPS || nal_type == MS_HEVC_NAL_PPS)
		array->array_completeness = ps_array_completeness;
	return 0;
}
static void hvcc_update_ptl(HEVCDecoderConfigurationRecord *hvcc, HVCCProfileTierLevel *ptl)
{

	hvcc->general_profile_space = ptl->profile_space;
	if (hvcc->general_tier_flag < ptl->tier_flag)
		hvcc->general_level_idc = ptl->level_idc;
	else
		hvcc->general_level_idc = FFMAX(hvcc->general_level_idc, ptl->level_idc);
	hvcc->general_tier_flag = FFMAX(hvcc->general_tier_flag, ptl->tier_flag);
	hvcc->general_profile_idc = FFMAX(hvcc->general_profile_idc, ptl->profile_idc);
	hvcc->general_profile_compatibility_flags &= ptl->profile_compatibility_flags;
	hvcc->general_constraint_indicator_flags &= ptl->constraint_indicator_flags;
}
static void hvcc_parse_ptl(GetBitContext *gb, HEVCDecoderConfigurationRecord *hvcc, unsigned int max_sub_layers_minus1)
{
	unsigned int i;
	HVCCProfileTierLevel general_ptl;
	uint8_t sub_layer_profile_present_flag[MS_HEVC_MAX_SUB_LAYERS];
	uint8_t sub_layer_level_present_flag[MS_HEVC_MAX_SUB_LAYERS];
	general_ptl.profile_space = get_bits(gb, 2);
	general_ptl.tier_flag = get_bits1(gb);
	general_ptl.profile_idc = get_bits(gb, 5);
	general_ptl.profile_compatibility_flags = get_bits_long(gb, 32);
	general_ptl.constraint_indicator_flags = get_bits64(gb, 48);
	general_ptl.level_idc = get_bits(gb, 8);
	for (i = 0; i < max_sub_layers_minus1; i++)
	{
		sub_layer_profile_present_flag[i] = get_bits1(gb);
		sub_layer_level_present_flag[i] = get_bits1(gb);
	}

	if (max_sub_layers_minus1 > 0)
		for (i = max_sub_layers_minus1; i < 8; i++)
			skip_bits(gb, 2);
	for (i = 0; i < max_sub_layers_minus1; i++)
	{
		if (sub_layer_profile_present_flag[i])
		{
			skip_bits_long(gb, 32);
			skip_bits_long(gb, 32);
			skip_bits(gb, 24);
		}
		if (sub_layer_level_present_flag[i])
			skip_bits(gb, 8);
	}
}
static void skip_sub_layer_hrd_parameters(GetBitContext *gb, unsigned int cpb_cnt_minus1, uint8_t sub_pic_hrd_params_present_flag)
{
	unsigned int i;
	for (i = 0; i <= cpb_cnt_minus1; i++)
	{
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);

		if (sub_pic_hrd_params_present_flag)
		{
			get_ue_golomb_long(gb);
			get_ue_golomb_long(gb);
		}

		skip_bits1(gb);
	}
}
static int skip_hrd_parameters(GetBitContext *gb, uint8_t cprms_present_flag, unsigned int max_sub_layers_minus1)
{
	unsigned int i;
	uint8_t sub_pic_hrd_params_present_flag = 0;
	uint8_t nal_hrd_parameters_present_flag = 0;
	uint8_t vcl_hrd_parameters_present_flag = 0;

	if (cprms_present_flag)
	{
		nal_hrd_parameters_present_flag = get_bits1(gb);
		vcl_hrd_parameters_present_flag = get_bits1(gb);

		if (nal_hrd_parameters_present_flag ||
			vcl_hrd_parameters_present_flag)
		{
			sub_pic_hrd_params_present_flag = get_bits1(gb);

			if (sub_pic_hrd_params_present_flag)
				skip_bits(gb, 19);
			skip_bits(gb, 8);

			if (sub_pic_hrd_params_present_flag)
				skip_bits(gb, 4);
			skip_bits(gb, 15);
		}
	}

	for (i = 0; i <= max_sub_layers_minus1; i++)
	{
		unsigned int cpb_cnt_minus1 = 0;
		uint8_t low_delay_hrd_flag = 0;
		uint8_t fixed_pic_rate_within_cvs_flag = 0;
		uint8_t fixed_pic_rate_general_flag = get_bits1(gb);

		if (!fixed_pic_rate_general_flag)
			fixed_pic_rate_within_cvs_flag = get_bits1(gb);

		if (fixed_pic_rate_within_cvs_flag)
			get_ue_golomb_long(gb);
		else
			low_delay_hrd_flag = get_bits1(gb);

		if (!low_delay_hrd_flag)
		{
			cpb_cnt_minus1 = get_ue_golomb_long(gb);
			if (cpb_cnt_minus1 > 31)
				return AVERROR_INVALIDDATA;
		}

		if (nal_hrd_parameters_present_flag)
			skip_sub_layer_hrd_parameters(gb, cpb_cnt_minus1,
				sub_pic_hrd_params_present_flag);

		if (vcl_hrd_parameters_present_flag)
			skip_sub_layer_hrd_parameters(gb, cpb_cnt_minus1,
				sub_pic_hrd_params_present_flag);
	}

	return 0;
}
static void skip_scaling_list_data(GetBitContext *gb)
{
	int i, j, k, num_coeffs;

	for (i = 0; i < 4; i++)
		for (j = 0; j < (i == 3 ? 2 : 6); j++)
			if (!get_bits1(gb))
				get_ue_golomb_long(gb);
			else
			{
				num_coeffs = FFMIN(64, 1 << (4 + (i << 1)));
				if (i > 1)
					get_se_golomb_long(gb);
				for (k = 0; k < num_coeffs; k++)
					get_se_golomb_long(gb);
			}
}
static void skip_timing_info(GetBitContext *gb)
{
	skip_bits_long(gb, 32);
	skip_bits_long(gb, 32);

	if (get_bits1(gb))
		get_ue_golomb_long(gb);
}
static int parse_rps(GetBitContext *gb, unsigned int rps_idx, unsigned int num_rps, unsigned int num_delta_pocs[MS_HEVC_MAX_SHORT_TERM_REF_PIC_SETS])
{
	unsigned int i;
	if (rps_idx && get_bits1(gb))
	{
		if (rps_idx >= num_rps)
			return AVERROR_INVALIDDATA;

		skip_bits1(gb);
		get_ue_golomb_long(gb);

		num_delta_pocs[rps_idx] = 0;
		for (i = 0; i <= num_delta_pocs[rps_idx - 1]; i++)
		{
			uint8_t use_delta_flag = 0;
			uint8_t used_by_curr_pic_flag = get_bits1(gb);
			if (!used_by_curr_pic_flag)
				use_delta_flag = get_bits1(gb);

			if (used_by_curr_pic_flag || use_delta_flag)
				num_delta_pocs[rps_idx]++;
		}
	}
	else
	{
		unsigned int num_negative_pics = get_ue_golomb_long(gb);
		unsigned int num_positive_pics = get_ue_golomb_long(gb);

		if ((num_positive_pics + (uint64_t)num_negative_pics) * 2 > get_bits_left(gb))
			return AVERROR_INVALIDDATA;

		num_delta_pocs[rps_idx] = num_negative_pics + num_positive_pics;

		for (i = 0; i < num_negative_pics; i++)
		{
			get_ue_golomb_long(gb);
			skip_bits1(gb);
		}

		for (i = 0; i < num_positive_pics; i++)
		{
			get_ue_golomb_long(gb);
			skip_bits1(gb);
		}
	}

	return 0;
}
static void hvcc_parse_vui(GetBitContext *gb, HEVCDecoderConfigurationRecord *hvcc, unsigned int max_sub_layers_minus1)
{
	unsigned int min_spatial_segmentation_idc;
	if (get_bits1(gb))
		if (get_bits(gb, 8) == 255)
			skip_bits_long(gb, 32);
	if (get_bits1(gb))
		skip_bits1(gb);
	if (get_bits1(gb))
	{
		skip_bits(gb, 4);
		if (get_bits1(gb))
			skip_bits(gb, 24);
	}
	if (get_bits1(gb))
	{
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
	}
	skip_bits(gb, 3);
	if (get_bits1(gb))
	{
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
	}
	if (get_bits1(gb))
	{
		skip_timing_info(gb);

		if (get_bits1(gb))
			skip_hrd_parameters(gb, 1, max_sub_layers_minus1);
	}
	if (get_bits1(gb))
	{
		skip_bits(gb, 3);
		min_spatial_segmentation_idc = get_ue_golomb_long(gb);
		hvcc->min_spatial_segmentation_idc = FFMIN(hvcc->min_spatial_segmentation_idc, min_spatial_segmentation_idc);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
	}
}
static int hvcc_parse_vps(GetBitContext *gb, HEVCDecoderConfigurationRecord *hvcc)
{
	unsigned int vps_max_sub_layers_minus1;
	skip_bits(gb, 12);
	vps_max_sub_layers_minus1 = get_bits(gb, 3);
	hvcc->numTemporalLayers = FFMAX(hvcc->numTemporalLayers, vps_max_sub_layers_minus1 + 1);
	skip_bits(gb, 17);
	hvcc_parse_ptl(gb, hvcc, vps_max_sub_layers_minus1);
	return 0;
}
static void skip_sub_layer_ordering_info(GetBitContext *gb)
{
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
}
static int hvcc_parse_sps(GetBitContext *gb, HEVCDecoderConfigurationRecord *hvcc)
{
	unsigned int i, sps_max_sub_layers_minus1, log2_max_pic_order_cnt_lsb_minus4;
	unsigned int num_short_term_ref_pic_sets, num_delta_pocs[MS_HEVC_MAX_SHORT_TERM_REF_PIC_SETS];
	skip_bits(gb, 4);
	sps_max_sub_layers_minus1 = get_bits(gb, 3);
	hvcc->numTemporalLayers = FFMAX(hvcc->numTemporalLayers, sps_max_sub_layers_minus1 + 1);
	hvcc->temporalIdNested = get_bits1(gb);
	hvcc_parse_ptl(gb, hvcc, sps_max_sub_layers_minus1);
	get_ue_golomb_long(gb);
	hvcc->chromaFormat = get_ue_golomb_long(gb);
	if (hvcc->chromaFormat == 3)
		skip_bits1(gb);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);

	if (get_bits1(gb))
	{
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
	}
	hvcc->bitDepthLumaMinus8 = get_ue_golomb_long(gb);
	hvcc->bitDepthChromaMinus8 = get_ue_golomb_long(gb);
	log2_max_pic_order_cnt_lsb_minus4 = get_ue_golomb_long(gb);

	i = get_bits1(gb) ? 0 : sps_max_sub_layers_minus1;
	for (; i <= sps_max_sub_layers_minus1; i++)
		skip_sub_layer_ordering_info(gb);

	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);

	if (get_bits1(gb) && get_bits1(gb))
		skip_scaling_list_data(gb);

	skip_bits1(gb);
	skip_bits1(gb);

	if (get_bits1(gb))
	{
		skip_bits(gb, 4);
		skip_bits(gb, 4);
		get_ue_golomb_long(gb);
		get_ue_golomb_long(gb);
		skip_bits1(gb);
	}

	num_short_term_ref_pic_sets = get_ue_golomb_long(gb);
	if (num_short_term_ref_pic_sets > MS_HEVC_MAX_SHORT_TERM_REF_PIC_SETS)
		return AVERROR_INVALIDDATA;

	for (i = 0; i < num_short_term_ref_pic_sets; i++)
	{
		int ret = parse_rps(gb, i, num_short_term_ref_pic_sets, num_delta_pocs);
		if (ret < 0)
			return ret;
	}

	if (get_bits1(gb))
	{
		unsigned num_long_term_ref_pics_sps = get_ue_golomb_long(gb);
		if (num_long_term_ref_pics_sps > 31U)
			return AVERROR_INVALIDDATA;
		for (i = 0; i < num_long_term_ref_pics_sps; i++)
		{
			int len = FFMIN(log2_max_pic_order_cnt_lsb_minus4 + 4, 16);
			skip_bits(gb, len);
			skip_bits1(gb);
		}
	}
	skip_bits1(gb);
	skip_bits1(gb);
	if (get_bits1(gb))
		hvcc_parse_vui(gb, hvcc, sps_max_sub_layers_minus1);
	return 0;
}
static int hvcc_parse_pps(GetBitContext *gb, HEVCDecoderConfigurationRecord *hvcc)
{
	uint8_t tiles_enabled_flag, entropy_coding_sync_enabled_flag;
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	skip_bits(gb, 7);
	get_ue_golomb_long(gb);
	get_ue_golomb_long(gb);
	get_se_golomb_long(gb);
	skip_bits(gb, 2);
	if (get_bits1(gb))
		get_ue_golomb_long(gb);
	get_se_golomb_long(gb);
	get_se_golomb_long(gb);
	skip_bits(gb, 4);
	tiles_enabled_flag = get_bits1(gb);
	entropy_coding_sync_enabled_flag = get_bits1(gb);
	if (entropy_coding_sync_enabled_flag && tiles_enabled_flag)
		hvcc->parallelismType = 0;
	else if (entropy_coding_sync_enabled_flag)
		hvcc->parallelismType = 3;
	else if (tiles_enabled_flag)
		hvcc->parallelismType = 2;
	else
		hvcc->parallelismType = 1;
	return 0;
}
static uint8_t *nal_unit_extract_rbsp(const uint8_t *src, uint32_t src_len, uint32_t *dst_len)
{
	uint8_t *dst;
	uint32_t i, len;

	dst = (uint8_t*)av_malloc(src_len + AV_INPUT_BUFFER_PADDING_SIZE);
	if (!dst)
		return NULL;

	/* NAL unit header (2 bytes) */
	i = len = 0;
	while (i < 2 && i < src_len)
		dst[len++] = src[i++];

	while (i + 2 < src_len)
		if (!src[i] && !src[i + 1] && src[i + 2] == 3) {
			dst[len++] = src[i++];
			dst[len++] = src[i++];
			i++; // remove emulation_prevention_three_byte
		}
		else
			dst[len++] = src[i++];

	while (i < src_len)
		dst[len++] = src[i++];

	memset(dst + len, 0, AV_INPUT_BUFFER_PADDING_SIZE);

	*dst_len = len;
	return dst;
}
static void nal_unit_parse_header(GetBitContext *gb, uint8_t *nal_type)
{
	skip_bits1(gb);
	*nal_type = get_bits(gb, 6);
	skip_bits(gb, 9);
}
static int hvcc_add_nal_unit(uint8_t *nal_buf, uint32_t nal_size, int ps_array_completeness, HEVCDecoderConfigurationRecord *hvcc)
{
	int ret = 0;
	GetBitContext gbc;
	uint8_t nal_type;
	uint8_t *rbsp_buf;
	uint32_t rbsp_size;

	rbsp_buf = nal_unit_extract_rbsp(nal_buf, nal_size, &rbsp_size);
	if (!rbsp_buf)
	{
		ret = AVERROR(ENOMEM);
		goto end;
	}

	ret = init_get_bits8(&gbc, rbsp_buf, rbsp_size);
	if (ret < 0)
		goto end;

	nal_unit_parse_header(&gbc, &nal_type);

	switch (nal_type)
	{
	case MS_HEVC_NAL_VPS:
	case MS_HEVC_NAL_SPS:
	case MS_HEVC_NAL_PPS:
	case MS_HEVC_NAL_SEI_PREFIX:
	case MS_HEVC_NAL_SEI_SUFFIX:
		ret = hvcc_array_add_nal_unit(nal_buf, nal_size, nal_type,
			ps_array_completeness, hvcc);
		if (ret < 0)
			goto end;
		else if (nal_type == MS_HEVC_NAL_VPS)
			ret = hvcc_parse_vps(&gbc, hvcc);
		else if (nal_type == MS_HEVC_NAL_SPS)
			ret = hvcc_parse_sps(&gbc, hvcc);
		else if (nal_type == MS_HEVC_NAL_PPS)
			ret = hvcc_parse_pps(&gbc, hvcc);
		if (ret < 0)
			goto end;
		break;
	default:
		ret = AVERROR_INVALIDDATA;
		goto end;
	}

end:
	av_free(rbsp_buf);
	return ret;
}
static int get_hvcc(struct HEVCDecoderConfigurationRecord* hvcc, uint8_t *data, size_t size)
{
	const uint8_t *nal_start, *nal_end;
	const uint8_t *end = data + size;
	int type;
	int ret = 0;
	nal_start = FindRange(data, end);
	while (true)
	{
		while (nal_start < end && !*(nal_start++));
		if (nal_start == end)
			break;
		GetBitContext gbc;
		nal_end = FindRange(nal_start, end);
		int size = nal_end - nal_start;
		ret = init_get_bits8(&gbc, nal_start, size);
		if (ret < 0)
			return ret;
		type = (nal_start[0] >> 1) & 0x3F;
		switch (type)
		{
		case MS_HEVC_NAL_VPS:
		case MS_HEVC_NAL_SPS:
		case MS_HEVC_NAL_PPS:
		case MS_HEVC_NAL_SEI_PREFIX:
		case MS_HEVC_NAL_SEI_SUFFIX:
			ret = hvcc_add_nal_unit((uint8_t*)nal_start, size, 1, hvcc);
			if (ret < 0)
				return ret;
			break;
		default:
			ret = AVERROR_INVALIDDATA;
			return ret;
		}
		nal_start = nal_end;
	}
	return ret;
}

//////////////////////////////////////////////////////////////////////////
typedef struct HEVCH2645NAL
{
	GetBitContext	gb;
	INT32			nal_unit_type;
	INT32			temporal_id;
	INT32			nuh_layer_id;
	INT32			ref_idc;
	INT32			skipped_bytes;
	INT32			skipped_bytes_pos_size;
	INT32* skipped_bytes_pos;
} HEVCH2645NAL;

static INT32 hevc_parse_nal_header(HEVCH2645NAL* nal)
{
	GetBitContext* gb = &nal->gb;
	if (get_bits1(gb) != 0)
		return AVERROR_INVALIDDATA;
	nal->nal_unit_type = get_bits(gb, 6);
	nal->nuh_layer_id = get_bits(gb, 6);
	nal->temporal_id = get_bits(gb, 3) - 1;
	if (nal->temporal_id < 0)
		return AVERROR_INVALIDDATA;
	return 0;
}

static const INT32 CodeSize(const UINT8* s, INT32 size)
{
	if (size < 4)
		return -1;
	if (s[0] != 0 || s[1] != 0)
		return -1;
	if (s[2] == 1)
		return 3;
	if ((s[2] == 0 && s[3] == 1))
		return 4;
	return -1;
}

static int decode_profile_tier_level(GetBitContext* gb, PTLCommon* ptl)
{
	int i;
	if (get_bits_left(gb) < (2 + 1 + 5 + 32 + 4 + 43 + 1))
		return -1;
	ptl->profile_space = get_bits(gb, 2);
	ptl->tier_flag = get_bits1(gb);
	ptl->profile_idc = get_bits(gb, 5);
	for (i = 0; i < 32; i++)
	{
		ptl->profile_compatibility_flag[i] = get_bits1(gb);

		if (ptl->profile_idc == 0 && i > 0 && ptl->profile_compatibility_flag[i])
			ptl->profile_idc = i;
	}
	ptl->progressive_source_flag = get_bits1(gb);
	ptl->interlaced_source_flag = get_bits1(gb);
	ptl->non_packed_constraint_flag = get_bits1(gb);
	ptl->frame_only_constraint_flag = get_bits1(gb);

#define check_profile_idc(idc) \
        ptl->profile_idc == idc || ptl->profile_compatibility_flag[idc]

	if (check_profile_idc(4) || check_profile_idc(5) || check_profile_idc(6) ||
		check_profile_idc(7) || check_profile_idc(8) || check_profile_idc(9) ||
		check_profile_idc(10)) {

		ptl->max_12bit_constraint_flag = get_bits1(gb);
		ptl->max_10bit_constraint_flag = get_bits1(gb);
		ptl->max_8bit_constraint_flag = get_bits1(gb);
		ptl->max_422chroma_constraint_flag = get_bits1(gb);
		ptl->max_420chroma_constraint_flag = get_bits1(gb);
		ptl->max_monochrome_constraint_flag = get_bits1(gb);
		ptl->intra_constraint_flag = get_bits1(gb);
		ptl->one_picture_only_constraint_flag = get_bits1(gb);
		ptl->lower_bit_rate_constraint_flag = get_bits1(gb);

		if (check_profile_idc(5) || check_profile_idc(9) || check_profile_idc(10)) {
			ptl->max_14bit_constraint_flag = get_bits1(gb);
			skip_bits_long(gb, 33); // XXX_reserved_zero_33bits[0..32]
		}
		else {
			skip_bits_long(gb, 34); // XXX_reserved_zero_34bits[0..33]
		}
	}
	else if (check_profile_idc(2))
	{
		skip_bits(gb, 7);
		ptl->one_picture_only_constraint_flag = get_bits1(gb);
		skip_bits_long(gb, 35); // XXX_reserved_zero_35bits[0..34]
	}
	else
	{
		skip_bits_long(gb, 43); // XXX_reserved_zero_43bits[0..42]
	}
	if (check_profile_idc(1) || check_profile_idc(2) || check_profile_idc(3) || check_profile_idc(4) || check_profile_idc(5) || check_profile_idc(9))
		ptl->inbld_flag = get_bits1(gb);
	else
		skip_bits1(gb);
#undef check_profile_idc

	return 0;
}

static int parse_ptl(GetBitContext* gb, PTL* ptl, int max_num_sub_layers)
{
	int i;
	if (decode_profile_tier_level(gb, &ptl->general_ptl) < 0 || get_bits_left(gb) < (8 + (8 * 2 * (max_num_sub_layers - 1 > 0))))
	{
		return -1;
	}
	ptl->general_ptl.level_idc = get_bits(gb, 8);
	for (i = 0; i < max_num_sub_layers - 1; i++)
	{
		ptl->sub_layer_profile_present_flag[i] = get_bits1(gb);
		ptl->sub_layer_level_present_flag[i] = get_bits1(gb);
	}
	if (max_num_sub_layers - 1 > 0)
		for (i = max_num_sub_layers - 1; i < 8; i++)
			skip_bits(gb, 2);
	for (i = 0; i < max_num_sub_layers - 1; i++)
	{
		if (ptl->sub_layer_profile_present_flag[i] &&
			decode_profile_tier_level(gb, &ptl->sub_layer_ptl[i]) < 0)
		{
			return -1;
		}
		if (ptl->sub_layer_level_present_flag[i])
		{
			if (get_bits_left(gb) < 8)
				return -1;
			else
				ptl->sub_layer_ptl[i].level_idc = get_bits(gb, 8);
		}
	}
	return 0;
}
static void parse_sublayer_hrd(GetBitContext* gb, unsigned int nb_cpb, int subpic_params_present)
{
	for (UINT32 i = 0; i < nb_cpb; i++)
	{
		get_ue_golomb_long(gb); // bit_rate_value_minus1
		get_ue_golomb_long(gb); // cpb_size_value_minus1
		if (subpic_params_present)
		{
			get_ue_golomb_long(gb); // cpb_size_du_value_minus1
			get_ue_golomb_long(gb); // bit_rate_du_value_minus1
		}
		skip_bits1(gb); // cbr_flag
	}
}

static int decode_hrd(GetBitContext* gb, int common_inf_present, int max_sublayers)
{
	int nal_params_present = 0, vcl_params_present = 0;
	int subpic_params_present = 0;
	int i;
	if (common_inf_present)
	{
		nal_params_present = get_bits1(gb);
		vcl_params_present = get_bits1(gb);
		if (nal_params_present || vcl_params_present)
		{
			subpic_params_present = get_bits1(gb);
			if (subpic_params_present) {
				skip_bits(gb, 8); // tick_divisor_minus2
				skip_bits(gb, 5); // du_cpb_removal_delay_increment_length_minus1
				skip_bits(gb, 1); // sub_pic_cpb_params_in_pic_timing_sei_flag
				skip_bits(gb, 5); // dpb_output_delay_du_length_minus1
			}
			skip_bits(gb, 4); // bit_rate_scale
			skip_bits(gb, 4); // cpb_size_scale
			if (subpic_params_present)
				skip_bits(gb, 4);  // cpb_size_du_scale
			skip_bits(gb, 5); // initial_cpb_removal_delay_length_minus1
			skip_bits(gb, 5); // au_cpb_removal_delay_length_minus1
			skip_bits(gb, 5); // dpb_output_delay_length_minus1
		}
	}

	for (i = 0; i < max_sublayers; i++)
	{
		int low_delay = 0;
		unsigned int nb_cpb = 1;
		int fixed_rate = get_bits1(gb);
		if (!fixed_rate)
			fixed_rate = get_bits1(gb);
		if (fixed_rate)
			get_ue_golomb_long(gb);  // elemental_duration_in_tc_minus1
		else
			low_delay = get_bits1(gb);
		if (!low_delay)
		{
			nb_cpb = get_ue_golomb_long(gb) + 1;
			if (nb_cpb < 1 || nb_cpb > 32)
			{
				return AVERROR_INVALIDDATA;
			}
		}
		if (nal_params_present)
			parse_sublayer_hrd(gb, nb_cpb, subpic_params_present);
		if (vcl_params_present)
			parse_sublayer_hrd(gb, nb_cpb, subpic_params_present);
	}
	return 0;
}
static void remove_pps(HEVCParamSets* s, int id)
{
	if (s->pps_list[id] && s->pps == (const HEVCPPS*)s->pps_list[id]->data)
		s->pps = NULL;
	av_buffer_unref(&s->pps_list[id]);
}

static void remove_sps(HEVCParamSets* s, int id)
{
	int i;
	if (s->sps_list[id]) {
		if (s->sps == (const HEVCSPS*)s->sps_list[id]->data)
			s->sps = NULL;

		/* drop all PPS that depend on this SPS */
		for (i = 0; i < FF_ARRAY_ELEMS(s->pps_list); i++)
			if (s->pps_list[i] && ((HEVCPPS*)s->pps_list[i]->data)->sps_id == id)
				remove_pps(s, i);

		av_assert0(!(s->sps_list[id] && s->sps == (HEVCSPS*)s->sps_list[id]->data));
	}
	av_buffer_unref(&s->sps_list[id]);
}
static void remove_vps(HEVCParamSets* s, int id)
{
	int i;
	if (s->vps_list[id]) {
		if (s->vps == (const HEVCVPS*)s->vps_list[id]->data)
			s->vps = NULL;

		for (i = 0; i < FF_ARRAY_ELEMS(s->sps_list); i++)
			if (s->sps_list[i] && ((HEVCSPS*)s->sps_list[i]->data)->vps_id == id)
				remove_sps(s, i);
	}
	av_buffer_unref(&s->vps_list[id]);
}



static const uint8_t ff_hevc_diag_scan4x4_x[16] =
{
		0, 0, 1, 0,
		1, 2, 0, 1,
		2, 3, 1, 2,
		3, 2, 3, 3,
};

static const uint8_t ff_hevc_diag_scan4x4_y[16] =
{
	0, 1, 0, 2,
	1, 0, 3, 2,
	1, 0, 3, 2,
	1, 3, 2, 3,
};
static const uint8_t ff_hevc_diag_scan8x8_x[64] = {
	0, 0, 1, 0,
	1, 2, 0, 1,
	2, 3, 0, 1,
	2, 3, 4, 0,
	1, 2, 3, 4,
	5, 0, 1, 2,
	3, 4, 5, 6,
	0, 1, 2, 3,
	4, 5, 6, 7,
	1, 2, 3, 4,
	5, 6, 7, 2,
	3, 4, 5, 6,
	7, 3, 4, 5,
	6, 7, 4, 5,
	6, 7, 5, 6,
	7, 6, 7, 7,
};
static const uint8_t ff_hevc_diag_scan8x8_y[64] = {
	0, 1, 0, 2,
	1, 0, 3, 2,
	1, 0, 4, 3,
	2, 1, 0, 5,
	4, 3, 2, 1,
	0, 6, 5, 4,
	3, 2, 1, 0,
	7, 6, 5, 4,
	3, 2, 1, 0,
	7, 6, 5, 4,
	3, 2, 1, 7,
	6, 5, 4, 3,
	2, 7, 6, 5,
	4, 3, 7, 6,
	5, 4, 7, 6,
	5, 7, 6, 7,
};
static int map_pixel_format(HEVCSPS* sps)
{
	const AVPixFmtDescriptor* desc;
	switch (sps->bit_depth)
	{
	case 8:
		if (sps->chroma_format_idc == 0) sps->pix_fmt = AV_PIX_FMT_GRAY8;
		if (sps->chroma_format_idc == 1) sps->pix_fmt = AV_PIX_FMT_YUV420P;
		if (sps->chroma_format_idc == 2) sps->pix_fmt = AV_PIX_FMT_YUV422P;
		if (sps->chroma_format_idc == 3) sps->pix_fmt = AV_PIX_FMT_YUV444P;
		break;
	case 9:
		if (sps->chroma_format_idc == 0) sps->pix_fmt = AV_PIX_FMT_GRAY9;
		if (sps->chroma_format_idc == 1) sps->pix_fmt = AV_PIX_FMT_YUV420P9;
		if (sps->chroma_format_idc == 2) sps->pix_fmt = AV_PIX_FMT_YUV422P9;
		if (sps->chroma_format_idc == 3) sps->pix_fmt = AV_PIX_FMT_YUV444P9;
		break;
	case 10:
		if (sps->chroma_format_idc == 0) sps->pix_fmt = AV_PIX_FMT_GRAY10;
		if (sps->chroma_format_idc == 1) sps->pix_fmt = AV_PIX_FMT_YUV420P10;
		if (sps->chroma_format_idc == 2) sps->pix_fmt = AV_PIX_FMT_YUV422P10;
		if (sps->chroma_format_idc == 3) sps->pix_fmt = AV_PIX_FMT_YUV444P10;
		break;
	case 12:
		if (sps->chroma_format_idc == 0) sps->pix_fmt = AV_PIX_FMT_GRAY12;
		if (sps->chroma_format_idc == 1) sps->pix_fmt = AV_PIX_FMT_YUV420P12;
		if (sps->chroma_format_idc == 2) sps->pix_fmt = AV_PIX_FMT_YUV422P12;
		if (sps->chroma_format_idc == 3) sps->pix_fmt = AV_PIX_FMT_YUV444P12;
		break;
	default:
		return AVERROR_INVALIDDATA;
	}
	desc = av_pix_fmt_desc_get(sps->pix_fmt);
	if (!desc)
		return AVERROR(EINVAL);
	sps->hshift[0] = sps->vshift[0] = 0;
	sps->hshift[2] = sps->hshift[1] = desc->log2_chroma_w;
	sps->vshift[2] = sps->vshift[1] = desc->log2_chroma_h;
	sps->pixel_shift = sps->bit_depth > 8;
	return 0;
}
static const uint8_t default_scaling_list_intra[] =
{
	16, 16, 16, 16, 17, 18, 21, 24,
	16, 16, 16, 16, 17, 19, 22, 25,
	16, 16, 17, 18, 20, 22, 25, 29,
	16, 16, 18, 21, 24, 27, 31, 36,
	17, 17, 20, 24, 30, 35, 41, 47,
	18, 19, 22, 27, 35, 44, 54, 65,
	21, 22, 25, 31, 41, 54, 70, 88,
	24, 25, 29, 36, 47, 65, 88, 115
};
static const uint8_t default_scaling_list_inter[] =
{
	16, 16, 16, 16, 17, 18, 20, 24,
	16, 16, 16, 17, 18, 20, 24, 25,
	16, 16, 17, 18, 20, 24, 25, 28,
	16, 17, 18, 20, 24, 25, 28, 33,
	17, 18, 20, 24, 25, 28, 33, 41,
	18, 20, 24, 25, 28, 33, 41, 54,
	20, 24, 25, 28, 33, 41, 54, 71,
	24, 25, 28, 33, 41, 54, 71, 91
};
static const AVRational vui_sar[] =
{
{  0,   1 },
{  1,   1 },
{ 12,  11 },
{ 10,  11 },
{ 16,  11 },
{ 40,  33 },
{ 24,  11 },
{ 20,  11 },
{ 32,  11 },
{ 80,  33 },
{ 18,  11 },
{ 15,  11 },
{ 64,  33 },
{ 160, 99 },
{  4,   3 },
{  3,   2 },
{  2,   1 },
};
static void set_default_scaling_list_data(ScalingList* sl)
{
	int matrixId;

	for (matrixId = 0; matrixId < 6; matrixId++)
	{
		// 4x4 default is 16
		memset(sl->sl[0][matrixId], 16, 16);
		sl->sl_dc[0][matrixId] = 16; // default for 16x16
		sl->sl_dc[1][matrixId] = 16; // default for 32x32
	}
	memcpy(sl->sl[1][0], default_scaling_list_intra, 64);
	memcpy(sl->sl[1][1], default_scaling_list_intra, 64);
	memcpy(sl->sl[1][2], default_scaling_list_intra, 64);
	memcpy(sl->sl[1][3], default_scaling_list_inter, 64);
	memcpy(sl->sl[1][4], default_scaling_list_inter, 64);
	memcpy(sl->sl[1][5], default_scaling_list_inter, 64);
	memcpy(sl->sl[2][0], default_scaling_list_intra, 64);
	memcpy(sl->sl[2][1], default_scaling_list_intra, 64);
	memcpy(sl->sl[2][2], default_scaling_list_intra, 64);
	memcpy(sl->sl[2][3], default_scaling_list_inter, 64);
	memcpy(sl->sl[2][4], default_scaling_list_inter, 64);
	memcpy(sl->sl[2][5], default_scaling_list_inter, 64);
	memcpy(sl->sl[3][0], default_scaling_list_intra, 64);
	memcpy(sl->sl[3][1], default_scaling_list_intra, 64);
	memcpy(sl->sl[3][2], default_scaling_list_intra, 64);
	memcpy(sl->sl[3][3], default_scaling_list_inter, 64);
	memcpy(sl->sl[3][4], default_scaling_list_inter, 64);
	memcpy(sl->sl[3][5], default_scaling_list_inter, 64);
}
static int scaling_list_data(GetBitContext* gb, ScalingList* sl, HEVCSPS* sps)
{
	uint8_t scaling_list_pred_mode_flag;
	int32_t scaling_list_dc_coef[2][6];
	int size_id, matrix_id, pos = 0;
	int i = 0;
	for (size_id = 0; size_id < 4; size_id++)
	{
		for (matrix_id = 0; matrix_id < 6; matrix_id += ((size_id == 3) ? 3 : 1))
		{
			scaling_list_pred_mode_flag = get_bits1(gb);
			if (!scaling_list_pred_mode_flag)
			{
				unsigned int delta = get_ue_golomb_long(gb);
				if (delta)
				{
					delta *= (size_id == 3) ? 3 : 1;
					if (matrix_id < delta)
						return AVERROR_INVALIDDATA;
					memcpy(sl->sl[size_id][matrix_id],
						sl->sl[size_id][matrix_id - delta],
						size_id > 0 ? 64 : 16);
					if (size_id > 1)
						sl->sl_dc[size_id - 2][matrix_id] = sl->sl_dc[size_id - 2][matrix_id - delta];
				}
			}
			else
			{
				int next_coef, coef_num;
				int32_t scaling_list_delta_coef;
				next_coef = 8;
				coef_num = FFMIN(64, 1 << (4 + (size_id << 1)));
				if (size_id > 1)
				{
					scaling_list_dc_coef[size_id - 2][matrix_id] = get_se_golomb(gb) + 8;
					next_coef = scaling_list_dc_coef[size_id - 2][matrix_id];
					sl->sl_dc[size_id - 2][matrix_id] = next_coef;
				}
				for (i = 0; i < coef_num; i++)
				{
					if (size_id == 0)
						pos = 4 * ff_hevc_diag_scan4x4_y[i] +
						ff_hevc_diag_scan4x4_x[i];
					else
						pos = 8 * ff_hevc_diag_scan8x8_y[i] +
						ff_hevc_diag_scan8x8_x[i];

					scaling_list_delta_coef = get_se_golomb(gb);
					next_coef = (next_coef + 256U + scaling_list_delta_coef) % 256;
					sl->sl[size_id][matrix_id][pos] = next_coef;
				}
			}
		}

	}
	if (sps->chroma_format_idc == 3)
	{
		for (i = 0; i < 64; i++) {
			sl->sl[3][1][i] = sl->sl[2][1][i];
			sl->sl[3][2][i] = sl->sl[2][2][i];
			sl->sl[3][4][i] = sl->sl[2][4][i];
			sl->sl[3][5][i] = sl->sl[2][5][i];
		}
		sl->sl_dc[1][1] = sl->sl_dc[0][1];
		sl->sl_dc[1][2] = sl->sl_dc[0][2];
		sl->sl_dc[1][4] = sl->sl_dc[0][4];
		sl->sl_dc[1][5] = sl->sl_dc[0][5];
	}
	return 0;
}
static int hevc_parse_short_term_rps(GetBitContext* gb, ShortTermRPS* rps, const HEVCSPS* sps, int is_slice_header)
{
	uint8_t rps_predict = 0;
	int delta_poc;
	int k0 = 0;
	int k1 = 0;
	int k = 0;
	UINT32 i = 0;
	if (rps != sps->st_rps && sps->nb_st_rps)
		rps_predict = get_bits1(gb);
	if (rps_predict)
	{
		const ShortTermRPS* rps_ridx;
		int delta_rps;
		unsigned abs_delta_rps;
		uint8_t use_delta_flag = 0;
		uint8_t delta_rps_sign;
		if (is_slice_header)
		{
			unsigned int delta_idx = get_ue_golomb_long(gb) + 1;
			if (delta_idx > sps->nb_st_rps)
				return AVERROR_INVALIDDATA;
			rps_ridx = &sps->st_rps[sps->nb_st_rps - delta_idx];
			rps->rps_idx_num_delta_pocs = rps_ridx->num_delta_pocs;
		}
		else
			rps_ridx = &sps->st_rps[rps - sps->st_rps - 1];

		delta_rps_sign = get_bits1(gb);
		abs_delta_rps = get_ue_golomb_long(gb) + 1;
		if (abs_delta_rps < 1 || abs_delta_rps > 32768)
			return AVERROR_INVALIDDATA;
		delta_rps = (1 - (delta_rps_sign << 1)) * abs_delta_rps;
		for (i = 0; i <= rps_ridx->num_delta_pocs; i++)
		{
			int used = rps->used[k] = get_bits1(gb);
			if (!used)
				use_delta_flag = get_bits1(gb);
			if (used || use_delta_flag)
			{
				if (i < rps_ridx->num_delta_pocs)
					delta_poc = delta_rps + rps_ridx->delta_poc[i];
				else
					delta_poc = delta_rps;
				rps->delta_poc[k] = delta_poc;
				if (delta_poc < 0)
					k0++;
				else
					k1++;
				k++;
			}
		}

		if (k >= FF_ARRAY_ELEMS(rps->used))
			return AVERROR_INVALIDDATA;
		rps->num_delta_pocs = k;
		rps->num_negative_pics = k0;
		if (rps->num_delta_pocs != 0)
		{
			int used, tmp;
			for (i = 1; i < rps->num_delta_pocs; i++)
			{
				delta_poc = rps->delta_poc[i];
				used = rps->used[i];
				for (k = i - 1; k >= 0; k--)
				{
					tmp = rps->delta_poc[k];
					if (delta_poc < tmp)
					{
						rps->delta_poc[k + 1] = tmp;
						rps->used[k + 1] = rps->used[k];
						rps->delta_poc[k] = delta_poc;
						rps->used[k] = used;
					}
				}
			}
		}
		if ((rps->num_negative_pics >> 1) != 0)
		{
			int used;
			k = rps->num_negative_pics - 1;
			for (i = 0; i < rps->num_negative_pics >> 1; i++)
			{
				delta_poc = rps->delta_poc[i];
				used = rps->used[i];
				rps->delta_poc[i] = rps->delta_poc[k];
				rps->used[i] = rps->used[k];
				rps->delta_poc[k] = delta_poc;
				rps->used[k] = used;
				k--;
			}
		}
	}
	else
	{
		unsigned int prev, nb_positive_pics;
		rps->num_negative_pics = get_ue_golomb_long(gb);
		nb_positive_pics = get_ue_golomb_long(gb);
		if (rps->num_negative_pics >= HEVC_MAX_REFS || nb_positive_pics >= HEVC_MAX_REFS)
			return AVERROR_INVALIDDATA;
		rps->num_delta_pocs = rps->num_negative_pics + nb_positive_pics;
		if (rps->num_delta_pocs)
		{
			prev = 0;
			for (i = 0; i < rps->num_negative_pics; i++)
			{
				delta_poc = get_ue_golomb_long(gb) + 1;
				if (delta_poc < 1 || delta_poc > 32768)
					return AVERROR_INVALIDDATA;
				prev -= delta_poc;
				rps->delta_poc[i] = prev;
				rps->used[i] = get_bits1(gb);
			}
			prev = 0;
			for (i = 0; i < nb_positive_pics; i++)
			{
				delta_poc = get_ue_golomb_long(gb) + 1;
				if (delta_poc < 1 || delta_poc > 32768)
					return AVERROR_INVALIDDATA;
				prev += delta_poc;
				rps->delta_poc[rps->num_negative_pics + i] = prev;
				rps->used[rps->num_negative_pics + i] = get_bits1(gb);
			}
		}
	}
	return 0;
}
static void parse_vui(GetBitContext* gb, int apply_defdispwin, HEVCSPS* sps)
{
	VUI backup_vui;
	VUI* vui = &sps->vui;
	GetBitContext backup;
	int sar_present, alt = 0;
	sar_present = get_bits1(gb);
	if (sar_present)
	{
		uint8_t sar_idx = get_bits(gb, 8);
		if (sar_idx < FF_ARRAY_ELEMS(vui_sar))
			vui->sar = vui_sar[sar_idx];
		else if (sar_idx == 255)
		{
			vui->sar.num = get_bits(gb, 16);
			vui->sar.den = get_bits(gb, 16);
		}
	}

	vui->overscan_info_present_flag = get_bits1(gb);
	if (vui->overscan_info_present_flag)
		vui->overscan_appropriate_flag = get_bits1(gb);
	vui->video_signal_type_present_flag = get_bits1(gb);
	if (vui->video_signal_type_present_flag)
	{
		vui->video_format = get_bits(gb, 3);
		vui->video_full_range_flag = get_bits1(gb);
		vui->colour_description_present_flag = get_bits1(gb);
		if (vui->video_full_range_flag && sps->pix_fmt == AV_PIX_FMT_YUV420P)
			sps->pix_fmt = AV_PIX_FMT_YUVJ420P;
		if (vui->colour_description_present_flag)
		{
			vui->colour_primaries = get_bits(gb, 8);
			vui->transfer_characteristic = get_bits(gb, 8);
			vui->matrix_coeffs = get_bits(gb, 8);
			if (!av_color_primaries_name((enum AVColorPrimaries)vui->colour_primaries))
				vui->colour_primaries = AVCOL_PRI_UNSPECIFIED;
			if (!av_color_transfer_name((enum AVColorTransferCharacteristic)vui->transfer_characteristic))
				vui->transfer_characteristic = AVCOL_TRC_UNSPECIFIED;
			if (!av_color_space_name((enum AVColorSpace)vui->matrix_coeffs))
				vui->matrix_coeffs = AVCOL_SPC_UNSPECIFIED;
			if (vui->matrix_coeffs == AVCOL_SPC_RGB)
			{
				switch (sps->pix_fmt)
				{
				case AV_PIX_FMT_YUV444P:
					sps->pix_fmt = AV_PIX_FMT_GBRP;
					break;
				case AV_PIX_FMT_YUV444P10:
					sps->pix_fmt = AV_PIX_FMT_GBRP10;
					break;
				case AV_PIX_FMT_YUV444P12:
					sps->pix_fmt = AV_PIX_FMT_GBRP12;
					break;
				}
			}
		}
	}
	vui->chroma_loc_info_present_flag = get_bits1(gb);
	if (vui->chroma_loc_info_present_flag)
	{
		vui->chroma_sample_loc_type_top_field = get_ue_golomb_long(gb);
		vui->chroma_sample_loc_type_bottom_field = get_ue_golomb_long(gb);
	}
	vui->neutra_chroma_indication_flag = get_bits1(gb);
	vui->field_seq_flag = get_bits1(gb);
	vui->frame_field_info_present_flag = get_bits1(gb);
	memcpy(&backup, gb, sizeof(backup));
	memcpy(&backup_vui, vui, sizeof(backup_vui));
	if (get_bits_left(gb) >= 68 && show_bits_long(gb, 21) == 0x100000)
		vui->default_display_window_flag = 0;
	else
		vui->default_display_window_flag = get_bits1(gb);
	if (vui->default_display_window_flag)
	{
		int vert_mult = hevc_sub_height_c[sps->chroma_format_idc];
		int horiz_mult = hevc_sub_width_c[sps->chroma_format_idc];
		vui->def_disp_win.left_offset = get_ue_golomb_long(gb) * horiz_mult;
		vui->def_disp_win.right_offset = get_ue_golomb_long(gb) * horiz_mult;
		vui->def_disp_win.top_offset = get_ue_golomb_long(gb) * vert_mult;
		vui->def_disp_win.bottom_offset = get_ue_golomb_long(gb) * vert_mult;
	}
timing_info:
	vui->vui_timing_info_present_flag = get_bits1(gb);
	if (vui->vui_timing_info_present_flag)
	{
		if (get_bits_left(gb) < 66 && !alt)
		{
			memcpy(vui, &backup_vui, sizeof(backup_vui));
			memcpy(gb, &backup, sizeof(backup));
			alt = 1;
			goto timing_info;
		}
		vui->vui_num_units_in_tick = get_bits_long(gb, 32);
		vui->vui_time_scale = get_bits_long(gb, 32);
		vui->vui_poc_proportional_to_timing_flag = get_bits1(gb);
		if (vui->vui_poc_proportional_to_timing_flag)
			vui->vui_num_ticks_poc_diff_one_minus1 = get_ue_golomb_long(gb);
		vui->vui_hrd_parameters_present_flag = get_bits1(gb);
		if (vui->vui_hrd_parameters_present_flag)
			decode_hrd(gb, 1, sps->max_sub_layers);
	}

	vui->bitstream_restriction_flag = get_bits1(gb);
	if (vui->bitstream_restriction_flag)
	{
		if (get_bits_left(gb) < 8 && !alt)
		{
			memcpy(vui, &backup_vui, sizeof(backup_vui));
			memcpy(gb, &backup, sizeof(backup));
			alt = 1;
			goto timing_info;
		}
		vui->tiles_fixed_structure_flag = get_bits1(gb);
		vui->motion_vectors_over_pic_boundaries_flag = get_bits1(gb);
		vui->restricted_ref_pic_lists_flag = get_bits1(gb);
		vui->min_spatial_segmentation_idc = get_ue_golomb_long(gb);
		vui->max_bytes_per_pic_denom = get_ue_golomb_long(gb);
		vui->max_bits_per_min_cu_denom = get_ue_golomb_long(gb);
		vui->log2_max_mv_length_horizontal = get_ue_golomb_long(gb);
		vui->log2_max_mv_length_vertical = get_ue_golomb_long(gb);
	}

	if (get_bits_left(gb) < 1 && !alt)
	{
		memcpy(vui, &backup_vui, sizeof(backup_vui));
		memcpy(gb, &backup, sizeof(backup));
		alt = 1;
		goto timing_info;
	}
}
static void hevc_pps_free(void* opaque, uint8_t* data)
{
	HEVCPPS* pps = (HEVCPPS*)data;
	av_freep(&pps->column_width);
	av_freep(&pps->row_height);
	av_freep(&pps->col_bd);
	av_freep(&pps->row_bd);
	av_freep(&pps->col_idxX);
	av_freep(&pps->ctb_addr_rs_to_ts);
	av_freep(&pps->ctb_addr_ts_to_rs);
	av_freep(&pps->tile_pos_rs);
	av_freep(&pps->tile_id);
	av_freep(&pps->min_tb_addr_zs_tab);
	av_freep(&pps);
}
static int pps_range_extensions(GetBitContext* gb, HEVCPPS* pps, HEVCSPS* sps)
{
	int i;
	if (pps->transform_skip_enabled_flag)
	{
		pps->log2_max_transform_skip_block_size = get_ue_golomb_long(gb) + 2;
	}
	pps->cross_component_prediction_enabled_flag = get_bits1(gb);
	pps->chroma_qp_offset_list_enabled_flag = get_bits1(gb);
	if (pps->chroma_qp_offset_list_enabled_flag)
	{
		pps->diff_cu_chroma_qp_offset_depth = get_ue_golomb_long(gb);
		pps->chroma_qp_offset_list_len_minus1 = get_ue_golomb_long(gb);
		if (pps->chroma_qp_offset_list_len_minus1 > 5)
			return AVERROR_INVALIDDATA;
		for (i = 0; i <= pps->chroma_qp_offset_list_len_minus1; i++)
		{
			pps->cb_qp_offset_list[i] = get_se_golomb_long(gb);
			pps->cr_qp_offset_list[i] = get_se_golomb_long(gb);
		}
	}
	pps->log2_sao_offset_scale_luma = get_ue_golomb_long(gb);
	pps->log2_sao_offset_scale_chroma = get_ue_golomb_long(gb);
	if (pps->log2_sao_offset_scale_luma > FFMAX(sps->bit_depth - 10, 0) || pps->log2_sao_offset_scale_chroma > FFMAX(sps->bit_depth_chroma - 10, 0))
		return AVERROR_INVALIDDATA;
	return(0);
}
static inline int setup_pps(GetBitContext* gb, HEVCPPS* pps, HEVCSPS* sps)
{
	int log2_diff;
	int pic_area_in_ctbs;
	int i, j, x, y, ctb_addr_rs, tile_id;
	pps->col_bd = (unsigned int*)av_malloc_array(pps->num_tile_columns + 1, sizeof(*pps->col_bd));
	pps->row_bd = (unsigned int*)av_malloc_array(pps->num_tile_rows + 1, sizeof(*pps->row_bd));
	pps->col_idxX = (int*)av_malloc_array(sps->ctb_width, sizeof(*pps->col_idxX));
	if (!pps->col_bd || !pps->row_bd || !pps->col_idxX)
		return AVERROR(ENOMEM);
	if (pps->uniform_spacing_flag)
	{
		if (!pps->column_width)
		{
			pps->column_width = (unsigned int*)av_malloc_array(pps->num_tile_columns, sizeof(*pps->column_width));
			pps->row_height = (unsigned int*)av_malloc_array(pps->num_tile_rows, sizeof(*pps->row_height));
		}
		if (!pps->column_width || !pps->row_height)
			return AVERROR(ENOMEM);
		for (i = 0; i < pps->num_tile_columns; i++)
			pps->column_width[i] = ((i + 1) * sps->ctb_width) / pps->num_tile_columns - (i * sps->ctb_width) / pps->num_tile_columns;
		for (i = 0; i < pps->num_tile_rows; i++)
			pps->row_height[i] = ((i + 1) * sps->ctb_height) / pps->num_tile_rows - (i * sps->ctb_height) / pps->num_tile_rows;
	}
	pps->col_bd[0] = 0;
	for (i = 0; i < pps->num_tile_columns; i++)
		pps->col_bd[i + 1] = pps->col_bd[i] + pps->column_width[i];
	pps->row_bd[0] = 0;
	for (i = 0; i < pps->num_tile_rows; i++)
		pps->row_bd[i + 1] = pps->row_bd[i] + pps->row_height[i];
	for (i = 0, j = 0; i < sps->ctb_width; i++)
	{
		if (i > pps->col_bd[j])
			j++;
		pps->col_idxX[i] = j;
	}
	pic_area_in_ctbs = sps->ctb_width * sps->ctb_height;
	pps->ctb_addr_rs_to_ts = (int*)av_malloc_array(pic_area_in_ctbs, sizeof(*pps->ctb_addr_rs_to_ts));
	pps->ctb_addr_ts_to_rs = (int*)av_malloc_array(pic_area_in_ctbs, sizeof(*pps->ctb_addr_ts_to_rs));
	pps->tile_id = (int*)av_malloc_array(pic_area_in_ctbs, sizeof(*pps->tile_id));
	pps->min_tb_addr_zs_tab = (int*)av_malloc_array((sps->tb_mask + 2) * (sps->tb_mask + 2), sizeof(*pps->min_tb_addr_zs_tab));
	if (!pps->ctb_addr_rs_to_ts || !pps->ctb_addr_ts_to_rs || !pps->tile_id || !pps->min_tb_addr_zs_tab)
		return AVERROR(ENOMEM);
	for (ctb_addr_rs = 0; ctb_addr_rs < pic_area_in_ctbs; ctb_addr_rs++)
	{
		int tb_x = ctb_addr_rs % sps->ctb_width;
		int tb_y = ctb_addr_rs / sps->ctb_width;
		int tile_x = 0;
		int tile_y = 0;
		int val = 0;
		for (i = 0; i < pps->num_tile_columns; i++)
		{
			if (tb_x < pps->col_bd[i + 1])
			{
				tile_x = i;
				break;
			}
		}

		for (i = 0; i < pps->num_tile_rows; i++)
		{
			if (tb_y < pps->row_bd[i + 1])
			{
				tile_y = i;
				break;
			}
		}
		for (i = 0; i < tile_x; i++)
			val += pps->row_height[tile_y] * pps->column_width[i];
		for (i = 0; i < tile_y; i++)
			val += sps->ctb_width * pps->row_height[i];
		val += (tb_y - pps->row_bd[tile_y]) * pps->column_width[tile_x] +
			tb_x - pps->col_bd[tile_x];
		pps->ctb_addr_rs_to_ts[ctb_addr_rs] = val;
		pps->ctb_addr_ts_to_rs[val] = ctb_addr_rs;
	}
	for (j = 0, tile_id = 0; j < pps->num_tile_rows; j++)
		for (i = 0; i < pps->num_tile_columns; i++, tile_id++)
			for (y = pps->row_bd[j]; y < pps->row_bd[j + 1]; y++)
				for (x = pps->col_bd[i]; x < pps->col_bd[i + 1]; x++)
					pps->tile_id[pps->ctb_addr_rs_to_ts[y * sps->ctb_width + x]] = tile_id;
	pps->tile_pos_rs = (int*)av_malloc_array(tile_id, sizeof(*pps->tile_pos_rs));
	if (!pps->tile_pos_rs)
		return AVERROR(ENOMEM);
	for (j = 0; j < pps->num_tile_rows; j++)
		for (i = 0; i < pps->num_tile_columns; i++)
			pps->tile_pos_rs[j * pps->num_tile_columns + i] =
			pps->row_bd[j] * sps->ctb_width + pps->col_bd[i];
	log2_diff = sps->log2_ctb_size - sps->log2_min_tb_size;
	pps->min_tb_addr_zs = &pps->min_tb_addr_zs_tab[1 * (sps->tb_mask + 2) + 1];
	for (y = 0; y < sps->tb_mask + 2; y++)
	{
		pps->min_tb_addr_zs_tab[y * (sps->tb_mask + 2)] = -1;
		pps->min_tb_addr_zs_tab[y] = -1;
	}
	for (y = 0; y < sps->tb_mask + 1; y++)
	{
		for (x = 0; x < sps->tb_mask + 1; x++)
		{
			int tb_x = x >> log2_diff;
			int tb_y = y >> log2_diff;
			int rs = sps->ctb_width * tb_y + tb_x;
			int val = pps->ctb_addr_rs_to_ts[rs] << (log2_diff * 2);
			for (i = 0; i < log2_diff; i++)
			{
				int m = 1 << i;
				val += (m & x ? m * m : 0) + (m & y ? 2 * m * m : 0);
			}
			pps->min_tb_addr_zs[y * (sps->tb_mask + 2) + x] = val;
		}
	}
	return 0;
}
static int hevc_parse_sps(HEVCSPS* sps, GetBitContext* gb, unsigned int* sps_id, int apply_defdispwin, AVBufferRef** vps_list)
{
	HEVCWindow* ow;
	int ret = 0;
	int log2_diff_max_min_transform_block_size;
	int bit_depth_chroma, start, vui_present, sublayer_ordering_info;
	int i;
	sps->vps_id = get_bits(gb, 4);
	if (vps_list && !vps_list[sps->vps_id])
		return AVERROR_INVALIDDATA;
	sps->max_sub_layers = get_bits(gb, 3) + 1;
	if (sps->max_sub_layers > HEVC_MAX_SUB_LAYERS)
		return AVERROR_INVALIDDATA;
	sps->temporal_id_nesting_flag = get_bits(gb, 1);
	if ((ret = parse_ptl(gb, &sps->ptl, sps->max_sub_layers)) < 0)
		return ret;
	*sps_id = get_ue_golomb_long(gb);
	if (*sps_id >= HEVC_MAX_SPS_COUNT)
		return AVERROR_INVALIDDATA;
	sps->chroma_format_idc = get_ue_golomb_long(gb);
	if (sps->chroma_format_idc > 3U)
		return AVERROR_INVALIDDATA;
	if (sps->chroma_format_idc == 3)
		sps->separate_colour_plane_flag = get_bits1(gb);
	if (sps->separate_colour_plane_flag)
		sps->chroma_format_idc = 0;
	sps->width = get_ue_golomb_long(gb);
	sps->height = get_ue_golomb_long(gb);
	if (get_bits1(gb))
	{
		int vert_mult = hevc_sub_height_c[sps->chroma_format_idc];
		int horiz_mult = hevc_sub_width_c[sps->chroma_format_idc];
		sps->pic_conf_win.left_offset = get_ue_golomb_long(gb) * horiz_mult;
		sps->pic_conf_win.right_offset = get_ue_golomb_long(gb) * horiz_mult;
		sps->pic_conf_win.top_offset = get_ue_golomb_long(gb) * vert_mult;
		sps->pic_conf_win.bottom_offset = get_ue_golomb_long(gb) * vert_mult;
		sps->output_window = sps->pic_conf_win;
	}
	sps->bit_depth = get_ue_golomb_long(gb) + 8;
	bit_depth_chroma = get_ue_golomb_long(gb) + 8;
	if (sps->chroma_format_idc && bit_depth_chroma != sps->bit_depth)
		return AVERROR_INVALIDDATA;
	sps->bit_depth_chroma = bit_depth_chroma;
	ret = map_pixel_format(sps);
	if (ret < 0)
		return ret;
	sps->log2_max_poc_lsb = get_ue_golomb_long(gb) + 4;
	if (sps->log2_max_poc_lsb > 16)
		return AVERROR_INVALIDDATA;
	sublayer_ordering_info = get_bits1(gb);
	start = sublayer_ordering_info ? 0 : sps->max_sub_layers - 1;
	for (i = start; i < sps->max_sub_layers; i++)
	{
		sps->temporal_layer[i].max_dec_pic_buffering = get_ue_golomb_long(gb) + 1;
		sps->temporal_layer[i].num_reorder_pics = get_ue_golomb_long(gb);
		sps->temporal_layer[i].max_latency_increase = get_ue_golomb_long(gb) - 1;
		if (sps->temporal_layer[i].max_dec_pic_buffering > (unsigned)HEVC_MAX_DPB_SIZE)
			return AVERROR_INVALIDDATA;
		if (sps->temporal_layer[i].num_reorder_pics > sps->temporal_layer[i].max_dec_pic_buffering - 1)
			sps->temporal_layer[i].max_dec_pic_buffering = sps->temporal_layer[i].num_reorder_pics + 1;
	}
	if (!sublayer_ordering_info)
	{
		for (int i = 0; i < start; i++)
		{
			sps->temporal_layer[i].max_dec_pic_buffering = sps->temporal_layer[start].max_dec_pic_buffering;
			sps->temporal_layer[i].num_reorder_pics = sps->temporal_layer[start].num_reorder_pics;
			sps->temporal_layer[i].max_latency_increase = sps->temporal_layer[start].max_latency_increase;
		}
	}

	sps->log2_min_cb_size = get_ue_golomb_long(gb) + 3;
	sps->log2_diff_max_min_coding_block_size = get_ue_golomb_long(gb);
	sps->log2_min_tb_size = get_ue_golomb_long(gb) + 2;
	log2_diff_max_min_transform_block_size = get_ue_golomb_long(gb);
	sps->log2_max_trafo_size = log2_diff_max_min_transform_block_size + sps->log2_min_tb_size;
	if (sps->log2_min_cb_size < 3 || sps->log2_min_cb_size > 30)
		return AVERROR_INVALIDDATA;
	if (sps->log2_diff_max_min_coding_block_size > 30)
		return AVERROR_INVALIDDATA;
	if (sps->log2_min_tb_size >= sps->log2_min_cb_size || sps->log2_min_tb_size < 2)
		return AVERROR_INVALIDDATA;
	if (log2_diff_max_min_transform_block_size < 0 || log2_diff_max_min_transform_block_size > 30)
		return AVERROR_INVALIDDATA;
	sps->max_transform_hierarchy_depth_inter = get_ue_golomb_long(gb);
	sps->max_transform_hierarchy_depth_intra = get_ue_golomb_long(gb);
	sps->scaling_list_enable_flag = get_bits1(gb);
	if (sps->scaling_list_enable_flag)
	{
		set_default_scaling_list_data(&sps->scaling_list);
		if (get_bits1(gb))
		{
			ret = scaling_list_data(gb, &sps->scaling_list, sps);
			if (ret < 0)
				return ret;
		}
	}
	sps->amp_enabled_flag = get_bits1(gb);
	sps->sao_enabled = get_bits1(gb);
	sps->pcm_enabled_flag = get_bits1(gb);
	if (sps->pcm_enabled_flag)
	{
		sps->pcm.bit_depth = get_bits(gb, 4) + 1;
		sps->pcm.bit_depth_chroma = get_bits(gb, 4) + 1;
		sps->pcm.log2_min_pcm_cb_size = get_ue_golomb_long(gb) + 3;
		sps->pcm.log2_max_pcm_cb_size = sps->pcm.log2_min_pcm_cb_size + get_ue_golomb_long(gb);
		if (FFMAX(sps->pcm.bit_depth, sps->pcm.bit_depth_chroma) > sps->bit_depth)
			return AVERROR_INVALIDDATA;
		sps->pcm.loop_filter_disable_flag = get_bits1(gb);
	}
	sps->nb_st_rps = get_ue_golomb_long(gb);
	if (sps->nb_st_rps > HEVC_MAX_SHORT_TERM_REF_PIC_SETS)
		return AVERROR_INVALIDDATA;
	for (i = 0; i < sps->nb_st_rps; i++)
	{
		if ((ret = hevc_parse_short_term_rps(gb, &sps->st_rps[i], sps, 0)) < 0)
			return ret;
	}

	sps->long_term_ref_pics_present_flag = get_bits1(gb);
	if (sps->long_term_ref_pics_present_flag)
	{
		sps->num_long_term_ref_pics_sps = get_ue_golomb_long(gb);
		if (sps->num_long_term_ref_pics_sps > HEVC_MAX_LONG_TERM_REF_PICS)
			return AVERROR_INVALIDDATA;
		for (i = 0; i < sps->num_long_term_ref_pics_sps; i++) {
			sps->lt_ref_pic_poc_lsb_sps[i] = get_bits(gb, sps->log2_max_poc_lsb);
			sps->used_by_curr_pic_lt_sps_flag[i] = get_bits1(gb);
		}
	}

	sps->sps_temporal_mvp_enabled_flag = get_bits1(gb);
	sps->sps_strong_intra_smoothing_enable_flag = get_bits1(gb);
	sps->vui.sar = AVRational{ 0, 1 };
	vui_present = get_bits1(gb);
	if (vui_present)
		parse_vui(gb, apply_defdispwin, sps);
	if (get_bits1(gb))
	{
		sps->sps_range_extension_flag = get_bits1(gb);
		skip_bits(gb, 7);
		if (sps->sps_range_extension_flag)
		{
			sps->transform_skip_rotation_enabled_flag = get_bits1(gb);
			sps->transform_skip_context_enabled_flag = get_bits1(gb);
			sps->implicit_rdpcm_enabled_flag = get_bits1(gb);
			sps->explicit_rdpcm_enabled_flag = get_bits1(gb);
			sps->extended_precision_processing_flag = get_bits1(gb);
			sps->intra_smoothing_disabled_flag = get_bits1(gb);
			sps->high_precision_offsets_enabled_flag = get_bits1(gb);
			sps->persistent_rice_adaptation_enabled_flag = get_bits1(gb);
			sps->cabac_bypass_alignment_enabled_flag = get_bits1(gb);
		}
	}
	if (apply_defdispwin)
	{
		sps->output_window.left_offset += sps->vui.def_disp_win.left_offset;
		sps->output_window.right_offset += sps->vui.def_disp_win.right_offset;
		sps->output_window.top_offset += sps->vui.def_disp_win.top_offset;
		sps->output_window.bottom_offset += sps->vui.def_disp_win.bottom_offset;
	}

	ow = &sps->output_window;
	if (ow->left_offset >= INT_MAX - ow->right_offset ||
		ow->top_offset >= INT_MAX - ow->bottom_offset ||
		ow->left_offset + ow->right_offset >= sps->width ||
		ow->top_offset + ow->bottom_offset >= sps->height)
	{
		memset(ow, 0, sizeof(*ow));
		memset(&sps->pic_conf_win, 0, sizeof(sps->pic_conf_win));
	}
	sps->log2_ctb_size = sps->log2_min_cb_size + sps->log2_diff_max_min_coding_block_size;
	sps->log2_min_pu_size = sps->log2_min_cb_size - 1;
	if (sps->log2_ctb_size > HEVC_MAX_LOG2_CTB_SIZE)
		return AVERROR_INVALIDDATA;
	if (sps->log2_ctb_size < 4)
		return AVERROR_INVALIDDATA;
	sps->ctb_width = (sps->width + (1 << sps->log2_ctb_size) - 1) >> sps->log2_ctb_size;
	sps->ctb_height = (sps->height + (1 << sps->log2_ctb_size) - 1) >> sps->log2_ctb_size;
	sps->ctb_size = sps->ctb_width * sps->ctb_height;
	sps->min_cb_width = sps->width >> sps->log2_min_cb_size;
	sps->min_cb_height = sps->height >> sps->log2_min_cb_size;
	sps->min_tb_width = sps->width >> sps->log2_min_tb_size;
	sps->min_tb_height = sps->height >> sps->log2_min_tb_size;
	sps->min_pu_width = sps->width >> sps->log2_min_pu_size;
	sps->min_pu_height = sps->height >> sps->log2_min_pu_size;
	sps->tb_mask = (1 << (sps->log2_ctb_size - sps->log2_min_tb_size)) - 1;
	sps->qp_bd_offset = 6 * (sps->bit_depth - 8);
	if (av_mod_uintp2(sps->width, sps->log2_min_cb_size) || av_mod_uintp2(sps->height, sps->log2_min_cb_size))
		return AVERROR_INVALIDDATA;
	if (sps->max_transform_hierarchy_depth_inter > static_cast<int>(sps->log2_ctb_size - sps->log2_min_tb_size))
		return AVERROR_INVALIDDATA;
	if (sps->max_transform_hierarchy_depth_intra > static_cast<int>(sps->log2_ctb_size - sps->log2_min_tb_size))
		return AVERROR_INVALIDDATA;
	if (sps->log2_max_trafo_size > FFMIN(sps->log2_ctb_size, 5))
		return AVERROR_INVALIDDATA;
	if (get_bits_left(gb) < 0)
		return AVERROR_INVALIDDATA;
	return 0;
}
static int hevc_parse_nal_vps(GetBitContext* gb, HEVCParamSets* ps)
{
	int i, j;
	int vps_id = 0;
	ptrdiff_t nal_size;
	HEVCVPS* vps;
	AVBufferRef* vps_buf = av_buffer_allocz(sizeof(*vps));
	if (!vps_buf)
		return AVERROR(ENOMEM);
	vps = (HEVCVPS*)vps_buf->data;
	nal_size = gb->buffer_end - gb->buffer;
	if (nal_size > sizeof(vps->data))
		vps->data_size = sizeof(vps->data);
	else
		vps->data_size = nal_size;
	memcpy(vps->data, gb->buffer, vps->data_size);
	vps_id = get_bits(gb, 4);
	if (get_bits(gb, 2) != 3)
		goto err;
	vps->vps_max_layers = get_bits(gb, 6) + 1;
	vps->vps_max_sub_layers = get_bits(gb, 3) + 1;
	vps->vps_temporal_id_nesting_flag = get_bits1(gb);
	if (get_bits(gb, 16) != 0xffff)
		goto err;
	if (vps->vps_max_sub_layers > HEVC_MAX_SUB_LAYERS)
		goto err;
	if (parse_ptl(gb, &vps->ptl, vps->vps_max_sub_layers) < 0)
		goto err;
	vps->vps_sub_layer_ordering_info_present_flag = get_bits1(gb);
	i = vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers - 1;
	for (; i < vps->vps_max_sub_layers; i++)
	{
		vps->vps_max_dec_pic_buffering[i] = get_ue_golomb_long(gb) + 1;
		vps->vps_num_reorder_pics[i] = get_ue_golomb_long(gb);
		vps->vps_max_latency_increase[i] = get_ue_golomb_long(gb) - 1;

		if (vps->vps_max_dec_pic_buffering[i] > HEVC_MAX_DPB_SIZE || !vps->vps_max_dec_pic_buffering[i])
			goto err;
	}
	vps->vps_max_layer_id = get_bits(gb, 6);
	vps->vps_num_layer_sets = get_ue_golomb_long(gb) + 1;
	if (vps->vps_num_layer_sets < 1 || vps->vps_num_layer_sets > 1024 || (vps->vps_num_layer_sets - 1LL) * (vps->vps_max_layer_id + 1LL) > get_bits_left(gb))
		goto err;
	for (i = 1; i < vps->vps_num_layer_sets; i++)
		for (j = 0; j <= vps->vps_max_layer_id; j++)
			skip_bits(gb, 1);
	vps->vps_timing_info_present_flag = get_bits1(gb);
	if (vps->vps_timing_info_present_flag)
	{
		vps->vps_num_units_in_tick = get_bits_long(gb, 32);
		vps->vps_time_scale = get_bits_long(gb, 32);
		vps->vps_poc_proportional_to_timing_flag = get_bits1(gb);
		if (vps->vps_poc_proportional_to_timing_flag)
			vps->vps_num_ticks_poc_diff_one = get_ue_golomb_long(gb) + 1;
		vps->vps_num_hrd_parameters = get_ue_golomb_long(gb);
		if (vps->vps_num_hrd_parameters > (unsigned)vps->vps_num_layer_sets)
			goto err;
		for (i = 0; i < vps->vps_num_hrd_parameters; i++)
		{
			int common_inf_present = 1;
			get_ue_golomb_long(gb);
			if (i)
				common_inf_present = get_bits1(gb);
			decode_hrd(gb, common_inf_present, vps->vps_max_sub_layers);
		}
	}
	get_bits1(gb);
	if (get_bits_left(gb) < 0)
	{
		if (ps->vps_list[vps_id])
			goto err;
	}
	if (ps->vps_list[vps_id] &&
		!memcmp(ps->vps_list[vps_id]->data, vps_buf->data, vps_buf->size))
	{
		av_buffer_unref(&vps_buf);
	}
	else
	{
		remove_vps(ps, vps_id);
		ps->vps_list[vps_id] = vps_buf;
	}
	return 0;
err:
	av_buffer_unref(&vps_buf);
	return AVERROR_INVALIDDATA;
}
static int hevc_parse_nal_sps(GetBitContext* gb, HEVCParamSets* ps, int apply_defdispwin)
{
	HEVCSPS* sps;
	AVBufferRef* sps_buf = av_buffer_allocz(sizeof(*sps));
	unsigned int sps_id;
	INT32 iRes = 0;
	ptrdiff_t nal_size;
	if (!sps_buf)
		return AVERROR(ENOMEM);
	sps = (HEVCSPS*)sps_buf->data;
	nal_size = gb->buffer_end - gb->buffer;
	if (nal_size > sizeof(sps->data))
	{
		sps->data_size = sizeof(sps->data);
	}
	else
	{
		sps->data_size = nal_size;
	}
	memcpy(sps->data, gb->buffer, sps->data_size);
	iRes = hevc_parse_sps(sps, gb, &sps_id, apply_defdispwin, ps->vps_list);
	if (iRes < 0)
	{
		av_buffer_unref(&sps_buf);
		return iRes;
	}
	if (ps->sps_list[sps_id] && !memcmp(ps->sps_list[sps_id]->data, sps_buf->data, sps_buf->size))
	{
		av_buffer_unref(&sps_buf);
	}
	else
	{
		remove_sps(ps, sps_id);
		ps->sps_list[sps_id] = sps_buf;
	}
	return 0;
}
static int hevc_parse_nal_pps(GetBitContext* gb, HEVCParamSets* ps)
{
	HEVCSPS* sps = NULL;
	int i, iRes = 0;
	unsigned int pps_id = 0;
	ptrdiff_t nal_size;
	unsigned log2_parallel_merge_level_minus2;
	AVBufferRef* pps_buf;
	HEVCPPS* pps = (HEVCPPS*)av_mallocz(sizeof(*pps));
	if (!pps)
		return AVERROR(ENOMEM);
	pps_buf = av_buffer_create((uint8_t*)pps, sizeof(*pps), hevc_pps_free, NULL, 0);
	if (!pps_buf)
	{
		av_freep(&pps);
		return AVERROR(ENOMEM);
	}
	nal_size = gb->buffer_end - gb->buffer;
	if (nal_size > sizeof(pps->data))
		pps->data_size = sizeof(pps->data);
	else
		pps->data_size = nal_size;
	memcpy(pps->data, gb->buffer, pps->data_size);
	pps->loop_filter_across_tiles_enabled_flag = 1;
	pps->num_tile_columns = 1;
	pps->num_tile_rows = 1;
	pps->uniform_spacing_flag = 1;
	pps->disable_dbf = 0;
	pps->beta_offset = 0;
	pps->tc_offset = 0;
	pps->log2_max_transform_skip_block_size = 2;
	pps_id = get_ue_golomb_long(gb);
	if (pps_id >= HEVC_MAX_PPS_COUNT)
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	pps->sps_id = get_ue_golomb_long(gb);
	if (pps->sps_id >= HEVC_MAX_SPS_COUNT)
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	if (!ps->sps_list[pps->sps_id])
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	sps = (HEVCSPS*)ps->sps_list[pps->sps_id]->data;
	pps->dependent_slice_segments_enabled_flag = get_bits1(gb);
	pps->output_flag_present_flag = get_bits1(gb);
	pps->num_extra_slice_header_bits = get_bits(gb, 3);
	pps->sign_data_hiding_flag = get_bits1(gb);
	pps->cabac_init_present_flag = get_bits1(gb);
	pps->num_ref_idx_l0_default_active = get_ue_golomb_long(gb) + 1;
	pps->num_ref_idx_l1_default_active = get_ue_golomb_long(gb) + 1;
	pps->pic_init_qp_minus26 = get_se_golomb(gb);
	pps->constrained_intra_pred_flag = get_bits1(gb);
	pps->transform_skip_enabled_flag = get_bits1(gb);
	pps->cu_qp_delta_enabled_flag = get_bits1(gb);
	pps->diff_cu_qp_delta_depth = 0;
	if (pps->cu_qp_delta_enabled_flag)
		pps->diff_cu_qp_delta_depth = get_ue_golomb_long(gb);
	if (pps->diff_cu_qp_delta_depth < 0 || (pps->diff_cu_qp_delta_depth > sps->log2_diff_max_min_coding_block_size))
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	pps->cb_qp_offset = get_se_golomb(gb);
	if (pps->cb_qp_offset < -12 || pps->cb_qp_offset > 12)
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	pps->cr_qp_offset = get_se_golomb(gb);
	if (pps->cr_qp_offset < -12 || pps->cr_qp_offset > 12)
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	pps->pic_slice_level_chroma_qp_offsets_present_flag = get_bits1(gb);
	pps->weighted_pred_flag = get_bits1(gb);
	pps->weighted_bipred_flag = get_bits1(gb);
	pps->transquant_bypass_enable_flag = get_bits1(gb);
	pps->tiles_enabled_flag = get_bits1(gb);
	pps->entropy_coding_sync_enabled_flag = get_bits1(gb);
	if (pps->tiles_enabled_flag)
	{
		int num_tile_columns_minus1 = get_ue_golomb(gb);
		int num_tile_rows_minus1 = get_ue_golomb(gb);
		if (num_tile_columns_minus1 < 0 || num_tile_columns_minus1 >= sps->ctb_width)
		{
			iRes = num_tile_columns_minus1 < 0 ? num_tile_columns_minus1 : AVERROR_INVALIDDATA;
			goto err;
		}
		if (num_tile_rows_minus1 < 0 || num_tile_rows_minus1 >= sps->ctb_height)
		{
			iRes = num_tile_rows_minus1 < 0 ? num_tile_rows_minus1 : AVERROR_INVALIDDATA;
			goto err;
		}
		pps->num_tile_columns = num_tile_columns_minus1 + 1;
		pps->num_tile_rows = num_tile_rows_minus1 + 1;
		pps->column_width = (unsigned int*)av_malloc_array(pps->num_tile_columns, sizeof(*pps->column_width));
		pps->row_height = (unsigned int*)av_malloc_array(pps->num_tile_rows, sizeof(*pps->row_height));
		if (!pps->column_width || !pps->row_height)
		{
			iRes = AVERROR(ENOMEM);
			goto err;
		}
		pps->uniform_spacing_flag = get_bits1(gb);
		if (!pps->uniform_spacing_flag)
		{
			uint64_t sum = 0;
			for (i = 0; i < pps->num_tile_columns - 1; i++)
			{
				pps->column_width[i] = get_ue_golomb_long(gb) + 1;
				sum += pps->column_width[i];
			}
			if (sum >= sps->ctb_width)
			{
				iRes = AVERROR_INVALIDDATA;
				goto err;
			}
			pps->column_width[pps->num_tile_columns - 1] = sps->ctb_width - sum;

			sum = 0;
			for (i = 0; i < pps->num_tile_rows - 1; i++)
			{
				pps->row_height[i] = get_ue_golomb_long(gb) + 1;
				sum += pps->row_height[i];
			}
			if (sum >= sps->ctb_height)
			{
				iRes = AVERROR_INVALIDDATA;
				goto err;
			}
			pps->row_height[pps->num_tile_rows - 1] = sps->ctb_height - sum;
		}
		pps->loop_filter_across_tiles_enabled_flag = get_bits1(gb);
	}
	pps->seq_loop_filter_across_slices_enabled_flag = get_bits1(gb);
	pps->deblocking_filter_control_present_flag = get_bits1(gb);
	if (pps->deblocking_filter_control_present_flag)
	{
		pps->deblocking_filter_override_enabled_flag = get_bits1(gb);
		pps->disable_dbf = get_bits1(gb);
		if (!pps->disable_dbf)
		{
			int beta_offset_div2 = get_se_golomb(gb);
			int tc_offset_div2 = get_se_golomb(gb);
			if (beta_offset_div2 < -6 || beta_offset_div2 > 6)
			{
				iRes = AVERROR_INVALIDDATA;
				goto err;
			}
			if (tc_offset_div2 < -6 || tc_offset_div2 > 6)
			{
				iRes = AVERROR_INVALIDDATA;
				goto err;
			}
			pps->beta_offset = 2 * beta_offset_div2;
			pps->tc_offset = 2 * tc_offset_div2;
		}
	}
	pps->scaling_list_data_present_flag = get_bits1(gb);
	if (pps->scaling_list_data_present_flag)
	{
		set_default_scaling_list_data(&pps->scaling_list);
		iRes = scaling_list_data(gb, &pps->scaling_list, sps);
		if (iRes < 0)
			goto err;
	}
	pps->lists_modification_present_flag = get_bits1(gb);
	log2_parallel_merge_level_minus2 = get_ue_golomb_long(gb);
	if (log2_parallel_merge_level_minus2 > sps->log2_ctb_size)
	{
		iRes = AVERROR_INVALIDDATA;
		goto err;
	}
	pps->log2_parallel_merge_level = log2_parallel_merge_level_minus2 + 2;
	pps->slice_header_extension_present_flag = get_bits1(gb);
	if (get_bits1(gb))
	{
		pps->pps_range_extensions_flag = get_bits1(gb);
		skip_bits(gb, 7);
		if (sps->ptl.general_ptl.profile_idc == FF_PROFILE_HEVC_REXT && pps->pps_range_extensions_flag)
		{
			if ((iRes = pps_range_extensions(gb, pps, sps)) < 0)
				goto err;
		}
	}
	iRes = setup_pps(gb, pps, sps);
	if (iRes < 0)
		goto err;
	if (get_bits_left(gb) < 0)
		goto err;
	remove_pps(ps, pps_id);
	ps->pps_list[pps_id] = pps_buf;
	return 0;
err:
	av_buffer_unref(&pps_buf);
	return iRes;
}
static int hevc_parse_slice_header(GetBitContext* gb, HEVCH2645NAL* s, HEVCParamSets* ps, HEVC_SLICE_HERDER* sh)
{
	INT32 i, iRes = 0;
	sh->first_slice_in_pic_flag = get_bits1(gb);
	sh->no_output_of_prior_pics_flag = 0;
	if (IS_IRAP(s))
		sh->no_output_of_prior_pics_flag = get_bits1(gb);
	sh->pps_id = get_ue_golomb_long(gb);
	if (sh->pps_id >= HEVC_MAX_PPS_COUNT || !ps->pps_list[sh->pps_id])
		return AVERROR_INVALIDDATA;
	if (!sh->first_slice_in_pic_flag && ps->pps != (HEVCPPS*)ps->pps_list[sh->pps_id]->data)
		return AVERROR_INVALIDDATA;
	ps->pps = (HEVCPPS*)ps->pps_list[sh->pps_id]->data;
	if (ps->sps != (HEVCSPS*)ps->sps_list[ps->pps->sps_id]->data)
	{
		const HEVCSPS* sps = (HEVCSPS*)ps->sps_list[ps->pps->sps_id]->data;
		const HEVCSPS* last_sps = ps->sps;
		if (last_sps && IS_IRAP(s) && s->nal_unit_type != MS_HEVC_NAL_CRA_NUT)
		{
			if (sps->width != last_sps->width || sps->height != last_sps->height ||
				sps->temporal_layer[sps->max_sub_layers - 1].max_dec_pic_buffering !=
				last_sps->temporal_layer[last_sps->max_sub_layers - 1].max_dec_pic_buffering)
				sh->no_output_of_prior_pics_flag = 0;
		}
	}
	sh->dependent_slice_segment_flag = 0;
	if (!sh->first_slice_in_pic_flag)
	{
		int slice_address_length;
		if (ps->pps->dependent_slice_segments_enabled_flag)
			sh->dependent_slice_segment_flag = get_bits1(gb);
		slice_address_length = av_ceil_log2(ps->sps->ctb_width * ps->sps->ctb_height);
		sh->slice_segment_addr = get_bitsz(gb, slice_address_length);
		if (sh->slice_segment_addr >= static_cast<UINT32>(ps->sps->ctb_width * ps->sps->ctb_height))
			return AVERROR_INVALIDDATA;
		if (!sh->dependent_slice_segment_flag)
		{
			sh->slice_addr = sh->slice_segment_addr;
		}
	}
	else
	{
		sh->slice_segment_addr = sh->slice_addr = 0;
	}

	if (!sh->dependent_slice_segment_flag)
	{
		for (i = 0; i < ps->pps->num_extra_slice_header_bits; i++)
			skip_bits(gb, 1);
		sh->slice_type = (enum HEVCSliceType)get_ue_golomb_long(gb);
		if (!(sh->slice_type == HEVC_SLICE_I ||
			sh->slice_type == HEVC_SLICE_P ||
			sh->slice_type == HEVC_SLICE_B))
			return AVERROR_INVALIDDATA;
		if (IS_IRAP(s) && sh->slice_type != HEVC_SLICE_I)
			return AVERROR_INVALIDDATA;
		sh->pic_output_flag = 1;
		if (ps->pps->output_flag_present_flag)
			sh->pic_output_flag = get_bits1(gb);
		if (ps->sps->separate_colour_plane_flag)
			sh->colour_plane_id = get_bits(gb, 2);
	}

	return 0;
}
static void hvcc_init(HEVCDecoderConfigurationRecord* hvcc)
{
	memset(hvcc, 0, sizeof(HEVCDecoderConfigurationRecord));
	hvcc->configurationVersion = 1;
	hvcc->lengthSizeMinusOne = 3;
	hvcc->general_profile_compatibility_flags = 0xffffffff;
	hvcc->general_constraint_indicator_flags = 0xffffffffffff;
	hvcc->min_spatial_segmentation_idc = MAX_SPATIAL_SEGMENTATION + 1;
}
static void hvcc_close(HEVCDecoderConfigurationRecord* hvcc)
{
	uint8_t i;
	for (i = 0; i < hvcc->numOfArrays; i++)
	{
		hvcc->array[i].numNalus = 0;
		av_freep(&hvcc->array[i].nalUnit);
		av_freep(&hvcc->array[i].nalUnitLength);
	}
	hvcc->numOfArrays = 0;
	av_freep(&hvcc->array);
}
#endif /*FF_HEVC_H*/