#ifndef FF_AVC_H
#define FF_AVC_H

#include <stdlib.h>
#include <windows.h>
extern "C"
{
#include <libavutil/avutil.h>
#include <libavutil/buffer.h>
#include <libavutil/avassert.h>
#include <libavutil/intreadwrite.h>
#include <libavcodec/avcodec.h>
}
#include "get_bits.h"
#include "common.h"
enum
{
	MS_H264_MAX_SPS_COUNT = 32,
	MS_H264_MAX_PPS_COUNT = 256,
	MS_H264_MAX_DPB_FRAMES = 16,
	MS_H264_MAX_REFS = 2 * MS_H264_MAX_DPB_FRAMES,
	MS_H264_MAX_RPLM_COUNT = MS_H264_MAX_REFS + 1,
	MS_H264_MAX_MMCO_COUNT = MS_H264_MAX_REFS * 2 + 3,
	MS_H264_MAX_SLICE_GROUPS = 8,
	MS_H264_MAX_CPB_CNT = 32,
	MS_H264_MAX_MB_PIC_SIZE = 139264,
	MS_H264_MAX_MB_WIDTH = 1055,
	MS_H264_MAX_MB_HEIGHT = 1055,
	MS_H264_MAX_WIDTH = MS_H264_MAX_MB_WIDTH * 16,
	MS_H264_MAX_HEIGHT = MS_H264_MAX_MB_HEIGHT * 16,
};

enum H264_SEI_TYPE
{
	H264_SEI_TYPE_BUFFERING_PERIOD = 0,
	H264_SEI_TYPE_PIC_TIMING = 1,
	H264_SEI_TYPE_PAN_SCAN_RECT = 2,
	H264_SEI_TYPE_FILLER_PAYLOAD = 3,
	H264_SEI_TYPE_USER_DATA_REGISTERED = 4,
	H264_SEI_TYPE_USER_DATA_UNREGISTERED = 5,
	H264_SEI_TYPE_RECOVERY_POINT = 6,
	H264_SEI_TYPE_FRAME_PACKING = 45,
	H264_SEI_TYPE_DISPLAY_ORIENTATION = 47,
	H264_SEI_TYPE_GREEN_METADATA = 56,
	H264_SEI_TYPE_MASTERING_DISPLAY_COLOUR_VOLUME = 137,
	H264_SEI_TYPE_ALTERNATIVE_TRANSFER = 147,

};

static const UINT8 UUID_ISO_IEC_11578[16] =
{
	  0xdc, 0x45, 0xe9, 0xbd, 0xe6, 0xd9, 0x48, 0xb7,
	  0x96, 0x2c, 0xd8, 0x20, 0xd9, 0x23, 0xee, 0xef
};
static inline UINT32 BITSTREAM_ENDIAN_FIX32(UINT32 x)
{
	return (x << 24) + ((x << 8) & 0xff0000) + ((x >> 8) & 0xff00) + (x >> 24);
}
static inline UINT64 BITSTREAM_ENDIAN_FIX64(UINT64 x)
{
	return BITSTREAM_ENDIAN_FIX32(x >> 32) + ((UINT64)BITSTREAM_ENDIAN_FIX32(x) << 32);
}


#define H264_MAX_PICTURE_COUNT	36
#define MAX_MMCO_COUNT			66
#define MAX_DELAYED_PIC_COUNT	16
#define QP_MAX_NUM				(51 + 6*6) 
#define MAX_SPS_COUNT			32
#define MAX_PPS_COUNT			256
#define MAX_LOG2_MAX_FRAME_NUM  (12 + 4)
#define MIN_LOG2_MAX_FRAME_NUM  4
#define EXTENDED_SAR			255

			/* picture type */
#define PICT_TOP_FIELD     1
#define PICT_BOTTOM_FIELD  2
#define PICT_FRAME         3

static const uint8_t default_scaling4[2][16] =
{
	{  6, 13, 20, 28, 13, 20, 28, 32,
	  20, 28, 32, 37, 28, 32, 37, 42 },
	{ 10, 14, 20, 24, 14, 20, 24, 27,
	  20, 24, 27, 30, 24, 27, 30, 34 }
};

static const uint8_t default_scaling8[2][64] =
{
	{  6, 10, 13, 16, 18, 23, 25, 27,
	  10, 11, 16, 18, 23, 25, 27, 29,
	  13, 16, 18, 23, 25, 27, 29, 31,
	  16, 18, 23, 25, 27, 29, 31, 33,
	  18, 23, 25, 27, 29, 31, 33, 36,
	  23, 25, 27, 29, 31, 33, 36, 38,
	  25, 27, 29, 31, 33, 36, 38, 40,
	  27, 29, 31, 33, 36, 38, 40, 42 },
	{  9, 13, 15, 17, 19, 21, 22, 24,
	  13, 13, 17, 19, 21, 22, 24, 25,
	  15, 17, 19, 21, 22, 24, 25, 27,
	  17, 19, 21, 22, 24, 25, 27, 28,
	  19, 21, 22, 24, 25, 27, 28, 30,
	  21, 22, 24, 25, 27, 28, 30, 32,
	  22, 24, 25, 27, 28, 30, 32, 33,
	  24, 25, 27, 28, 30, 32, 33, 35 }
};

/* maximum number of MBs in the DPB for a given level */
static const int level_max_dpb_mbs[][2] =
{
	{ 10, 396       },
	{ 11, 900       },
	{ 12, 2376      },
	{ 13, 2376      },
	{ 20, 2376      },
	{ 21, 4752      },
	{ 22, 8100      },
	{ 30, 8100      },
	{ 31, 18000     },
	{ 32, 20480     },
	{ 40, 32768     },
	{ 41, 32768     },
	{ 42, 34816     },
	{ 50, 110400    },
	{ 51, 184320    },
	{ 52, 184320    },
};
static const AVRational h264_pixel_aspect[17] =
{
	{   0,  1 },
	{   1,  1 },
	{  12, 11 },
	{  10, 11 },
	{  16, 11 },
	{  40, 33 },
	{  24, 11 },
	{  20, 11 },
	{  32, 11 },
	{  80, 33 },
	{  18, 11 },
	{  15, 11 },
	{  64, 33 },
	{ 160, 99 },
	{   4,  3 },
	{   3,  2 },
	{   2,  1 },
};

static const uint8_t zigzag_direct[64] =
{
	0,   1,  8, 16,  9,  2,  3, 10,
	17, 24, 32, 25, 18, 11,  4,  5,
	12, 19, 26, 33, 40, 48, 41, 34,
	27, 20, 13,  6,  7, 14, 21, 28,
	35, 42, 49, 56, 57, 50, 43, 36,
	29, 22, 15, 23, 30, 37, 44, 51,
	58, 59, 52, 45, 38, 31, 39, 46,
	53, 60, 61, 54, 47, 55, 62, 63
};

static const uint8_t zigzag_scan[16 + 1] =
{
	0 + 0 * 4, 1 + 0 * 4, 0 + 1 * 4, 0 + 2 * 4,
	1 + 1 * 4, 2 + 0 * 4, 3 + 0 * 4, 2 + 1 * 4,
	1 + 2 * 4, 0 + 3 * 4, 1 + 3 * 4, 2 + 2 * 4,
	3 + 1 * 4, 3 + 2 * 4, 2 + 3 * 4, 3 + 3 * 4,
};
static const uint8_t h264_dequant4_coeff_init[6][3] =
{
	{ 10, 13, 16 },
	{ 11, 14, 18 },
	{ 13, 16, 20 },
	{ 14, 18, 23 },
	{ 16, 20, 25 },
	{ 18, 23, 29 },
};

static  const uint8_t h264_dequant8_coeff_init_scan[16] =
{
	0, 3, 4, 3, 3, 1, 5, 1, 4, 5, 2, 5, 3, 1, 5, 1
};

static  const uint8_t h264_dequant8_coeff_init[6][6] =
{
	{ 20, 18, 32, 19, 25, 24 },
	{ 22, 19, 35, 21, 28, 26 },
	{ 26, 23, 42, 24, 33, 31 },
	{ 28, 25, 45, 26, 35, 33 },
	{ 32, 28, 51, 30, 40, 38 },
	{ 36, 32, 58, 34, 46, 43 },
};

static const uint8_t h264_quant_rem6[QP_MAX_NUM + 1] =
{
	0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2,
	3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5,
	0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2,
	3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5,
	0, 1, 2, 3,
};

static const uint8_t h264_quant_div6[QP_MAX_NUM + 1] =
{
	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3,  3,  3,
	3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,  6,  6,
	7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10,
   10,10,10,11,11,11,11,11,11,12,12,12,12,12,12,13,13,13, 13, 13, 13,
   14,14,14,14,
};

#define QP(qP, depth) ((qP) + 6 * ((depth) - 8))

#define CHROMA_QP_TABLE_END(d)                                          \
    QP(0,  d), QP(1,  d), QP(2,  d), QP(3,  d), QP(4,  d), QP(5,  d),   \
    QP(6,  d), QP(7,  d), QP(8,  d), QP(9,  d), QP(10, d), QP(11, d),   \
    QP(12, d), QP(13, d), QP(14, d), QP(15, d), QP(16, d), QP(17, d),   \
    QP(18, d), QP(19, d), QP(20, d), QP(21, d), QP(22, d), QP(23, d),   \
    QP(24, d), QP(25, d), QP(26, d), QP(27, d), QP(28, d), QP(29, d),   \
    QP(29, d), QP(30, d), QP(31, d), QP(32, d), QP(32, d), QP(33, d),   \
    QP(34, d), QP(34, d), QP(35, d), QP(35, d), QP(36, d), QP(36, d),   \
    QP(37, d), QP(37, d), QP(37, d), QP(38, d), QP(38, d), QP(38, d),   \
    QP(39, d), QP(39, d), QP(39, d), QP(39, d)

static const uint8_t h264_chroma_qp[7][QP_MAX_NUM + 1] = {
	{ CHROMA_QP_TABLE_END(8) },
	{ 0, 1, 2, 3, 4, 5,
	  CHROMA_QP_TABLE_END(9) },
	{ 0, 1, 2, 3,  4,  5,
	  6, 7, 8, 9, 10, 11,
	  CHROMA_QP_TABLE_END(10) },
	{ 0,  1, 2, 3,  4,  5,
	  6,  7, 8, 9, 10, 11,
	  12,13,14,15, 16, 17,
	  CHROMA_QP_TABLE_END(11) },
	{ 0,  1, 2, 3,  4,  5,
	  6,  7, 8, 9, 10, 11,
	  12,13,14,15, 16, 17,
	  18,19,20,21, 22, 23,
	  CHROMA_QP_TABLE_END(12) },
	{ 0,  1, 2, 3,  4,  5,
	  6,  7, 8, 9, 10, 11,
	  12,13,14,15, 16, 17,
	  18,19,20,21, 22, 23,
	  24,25,26,27, 28, 29,
	  CHROMA_QP_TABLE_END(13) },
	{ 0,  1, 2, 3,  4,  5,
	  6,  7, 8, 9, 10, 11,
	  12,13,14,15, 16, 17,
	  18,19,20,21, 22, 23,
	  24,25,26,27, 28, 29,
	  30,31,32,33, 34, 35,
	  CHROMA_QP_TABLE_END(14) },
};

static const UINT8 FF_PICTURE_TYPES[5] =
{
	AV_PICTURE_TYPE_P,
	AV_PICTURE_TYPE_B,
	AV_PICTURE_TYPE_I,
	AV_PICTURE_TYPE_SP,
	AV_PICTURE_TYPE_SI
};

typedef struct tagH264PredWeightTable
{
	int use_weight;
	int use_weight_chroma;
	int luma_log2_weight_denom;
	int chroma_log2_weight_denom;
	int luma_weight_flag[2];    ///< 7.4.3.2 luma_weight_lX_flag
	int chroma_weight_flag[2];  ///< 7.4.3.2 chroma_weight_lX_flag
	// The following 2 can be changed to int8_t but that causes a 10 CPU cycles speed loss
	int luma_weight[48][2][2];
	int chroma_weight[48][2][2][2];
	int implicit_weight[48][48][2];
} H264PredWeightTable;

typedef struct tagH264POCContext
{
	int poc_lsb;
	int poc_msb;
	int delta_poc_bottom;
	int delta_poc[2];
	int frame_num;
	int prev_poc_msb;           ///< poc_msb of the last reference pic for POC type 0
	int prev_poc_lsb;           ///< poc_lsb of the last reference pic for POC type 0
	int frame_num_offset;       ///< for POC type 2
	int prev_frame_num_offset;  ///< for POC type 2
	int prev_frame_num;         ///< frame_num of the last pic for POC type 1/2
} H264POCContext;

typedef struct tagSPS
{
	unsigned int sps_id;
	int profile_idc;
	int level_idc;
	int chroma_format_idc;
	int transform_bypass;              ///< qpprime_y_zero_transform_bypass_flag
	int log2_max_frame_num;            ///< log2_max_frame_num_minus4 + 4
	int poc_type;                      ///< pic_order_cnt_type
	int log2_max_poc_lsb;              ///< log2_max_pic_order_cnt_lsb_minus4
	int delta_pic_order_always_zero_flag;
	int offset_for_non_ref_pic;
	int offset_for_top_to_bottom_field;
	int poc_cycle_length;              ///< num_ref_frames_in_pic_order_cnt_cycle
	int ref_frame_count;               ///< num_ref_frames
	int gaps_in_frame_num_allowed_flag;
	int mb_width;                      ///< pic_width_in_mbs_minus1 + 1
	///< (pic_height_in_map_units_minus1 + 1) * (2 - frame_mbs_only_flag)
	int mb_height;
	int frame_mbs_only_flag;
	int mb_aff;                        ///< mb_adaptive_frame_field_flag
	int direct_8x8_inference_flag;
	int crop;                          ///< frame_cropping_flag
	/* those 4 are already in luma samples */
	unsigned int crop_left;            ///< frame_cropping_rect_left_offset
	unsigned int crop_right;           ///< frame_cropping_rect_right_offset
	unsigned int crop_top;             ///< frame_cropping_rect_top_offset
	unsigned int crop_bottom;          ///< frame_cropping_rect_bottom_offset
	int vui_parameters_present_flag;
	AVRational sar;
	int video_signal_type_present_flag;
	int full_range;
	int colour_description_present_flag;
	enum AVColorPrimaries color_primaries;
	enum AVColorTransferCharacteristic color_trc;
	enum AVColorSpace colorspace;
	enum AVChromaLocation chroma_location;
	int timing_info_present_flag;
	uint32_t num_units_in_tick;
	uint32_t time_scale;
	int fixed_frame_rate_flag;
	int32_t offset_for_ref_frame[256];
	int bitstream_restriction_flag;
	int num_reorder_frames;
	int scaling_matrix_present;
	uint8_t scaling_matrix4[6][16];
	uint8_t scaling_matrix8[6][64];
	int nal_hrd_parameters_present_flag;
	int vcl_hrd_parameters_present_flag;
	int pic_struct_present_flag;
	int time_offset_length;
	int cpb_cnt;                          ///< See H.264 E.1.2
	int initial_cpb_removal_delay_length; ///< initial_cpb_removal_delay_length_minus1 + 1
	int cpb_removal_delay_length;         ///< cpb_removal_delay_length_minus1 + 1
	int dpb_output_delay_length;          ///< dpb_output_delay_length_minus1 + 1
	int bit_depth_luma;                   ///< bit_depth_luma_minus8 + 8
	int bit_depth_chroma;                 ///< bit_depth_chroma_minus8 + 8
	int residual_color_transform_flag;    ///< residual_colour_transform_flag
	int constraint_set_flags;             ///< constraint_set[0-3]_flag
	uint8_t data[4096];
	size_t data_size;
}SPS;

typedef struct tagPPS
{
	unsigned int sps_id;
	int cabac;                  ///< entropy_coding_mode_flag
	int pic_order_present;      ///< pic_order_present_flag
	int slice_group_count;      ///< num_slice_groups_minus1 + 1
	int mb_slice_group_map_type;
	unsigned int ref_count[2];  ///< num_ref_idx_l0/1_active_minus1 + 1
	int weighted_pred;          ///< weighted_pred_flag
	int weighted_bipred_idc;
	int init_qp;                ///< pic_init_qp_minus26 + 26
	int init_qs;                ///< pic_init_qs_minus26 + 26
	int chroma_qp_index_offset[2];
	int deblocking_filter_parameters_present; ///< deblocking_filter_parameters_present_flag
	int constrained_intra_pred;     ///< constrained_intra_pred_flag
	int redundant_pic_cnt_present;  ///< redundant_pic_cnt_present_flag
	int transform_8x8_mode;         ///< transform_8x8_mode_flag
	uint8_t scaling_matrix4[6][16];
	uint8_t scaling_matrix8[6][64];
	uint8_t chroma_qp_table[2][QP_MAX_NUM + 1];  ///< pre-scaled (with chroma_qp_index_offset) version of qp_table
	int chroma_qp_diff;
	uint8_t data[4096];
	size_t data_size;
	uint32_t dequant4_buffer[6][QP_MAX_NUM + 1][16];
	uint32_t dequant8_buffer[6][QP_MAX_NUM + 1][64];
	uint32_t(*dequant4_coeff[6])[16];
	uint32_t(*dequant8_coeff[6])[64];
	AVBufferRef *sps_ref;
	const SPS   *sps;
}PPS;

typedef struct tagH264ParamSets
{
	AVBufferRef *sps_list[MAX_SPS_COUNT];
	AVBufferRef *pps_list[MAX_PPS_COUNT];
	AVBufferRef *pps_ref;
	/* currently active parameters sets */
	const PPS *pps;
	const SPS *sps;
	int overread_warning_printed[2];
}H264ParamSets;

typedef struct tagH2645NAL
{
	uint8_t *rbsp_buffer;
	int size;
	const uint8_t *data;
	/**
	 * Size, in bits, of just the data, excluding the stop bit and any trailing
	 * padding. I.e. what HEVC calls SODB.
	 */
	int size_bits;
	int raw_size;
	const uint8_t *raw_data;
	GetBitContext gb;
	/**
	 * NAL unit type
	 */
	int type;
	/**
	 * HEVC only, nuh_temporal_id_plus_1 - 1
	 */
	int temporal_id;
	/*
	 * HEVC only, identifier of layer to which nal unit belongs
	 */
	int nuh_layer_id;
	int skipped_bytes;
	int skipped_bytes_pos_size;
	int *skipped_bytes_pos;
	/**
	 * H.264 only, nal_ref_idc
	 */
	int ref_idc;
} H2645NAL;

typedef struct tagH2645RBSP
{
	uint8_t *rbsp_buffer;
	AVBufferRef *rbsp_buffer_ref;
	int rbsp_buffer_alloc_size;
	int rbsp_buffer_size;
} H2645RBSP;


typedef struct
{
uint8_t id;
uint8_t profile_idc;
uint8_t level_idc;
uint8_t constraint_set_flags;
uint8_t chroma_format_idc;
uint8_t bit_depth_luma;
uint8_t bit_depth_chroma;
uint8_t frame_mbs_only_flag;
AVRational sar;
} H264SPS;



static const BOOL  CheckCode(const UINT8* s, INT32 size)
{
	if (size < 4)
		return FALSE;
	if (s[0] != 0 || s[1] != 0)
		return FALSE;
	return s[2] == 1 || (s[2] == 0 && s[3] == 1);
}
FFFRAMEWORK_API const INT32 FindCodeSize(const UINT8* s, UINT32 size);
static const UINT8* FindCode(const UINT8* s, const UINT8* end)
{
	const UINT8* a = s + 4 - ((INT_PTR)s & 3);
	for (end -= 3; s < a && s < end; s++)
	{
		if (s[0] == 0 && s[1] == 0 && s[2] == 1)
			return s;
	}
	for (end -= 3; s < end; s += 4)
	{
		UINT32 x = *(const UINT32*)s;
		if ((x - 0x01010101) & (~x) & 0x80808080)
		{
			if (s[1] == 0)
			{
				if (s[0] == 0 && s[2] == 1)
					return s;
				if (s[2] == 0 && s[3] == 1)
					return s + 1;
			}
			if (s[3] == 0)
			{
				if (s[2] == 0 && s[4] == 1)
					return s + 2;
				if (s[4] == 0 && s[5] == 1)
					return s + 3;
			}
		}
	}
	for (end += 3; s < end; s++)
	{
		if (s[0] == 0 && s[1] == 0 && s[2] == 1)
			return s;
	}

	return end + 3;
}

FFFRAMEWORK_API const UINT8* FindRange(const UINT8* s, const UINT8* end);

static UINT8* FindCodeW(UINT8* s, UINT8* end)
{
	UINT8* a = s + 4 - ((INT_PTR)s & 3);
	for (end -= 3; s < a && s < end; s++)
	{
		if (s[0] == 0 && s[1] == 0 && s[2] == 1)
			return s;
	}
	for (end -= 3; s < end; s += 4)
	{
		UINT32 x = *(const UINT32*)s;
		if ((x - 0x01010101) & (~x) & 0x80808080)
		{
			if (s[1] == 0)
			{
				if (s[0] == 0 && s[2] == 1)
					return s;
				if (s[2] == 0 && s[3] == 1)
					return s + 1;
			}
			if (s[3] == 0)
			{
				if (s[2] == 0 && s[4] == 1)
					return s + 2;
				if (s[4] == 0 && s[5] == 1)
					return s + 3;
			}
		}
	}
	for (end += 3; s < end; s++)
	{
		if (s[0] == 0 && s[1] == 0 && s[2] == 1)
			return s;
	}

	return end + 3;
}
static UINT8* FindRangeW(UINT8* s, UINT8* end)
{
	UINT8* out = FindCodeW(s, end);
	if (s < out && out < end && !out[-1])
		out--;
	return out;
}
#endif