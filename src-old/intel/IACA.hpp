#ifndef IACA_HEADER_FOR_IMEN
#define IACA_HEADER_FOR_IMEN

    //#define IACA_MARKS_OFF
    #ifdef IACA_MARKS_OFF

        #define IACA_START
        #define IACA_END
        #define IACA_MSC64_START
        #define IACA_MSC64_END

    #else
        #if defined (__GNUC__)
            #define IACA_SSC_MARK( MARK_ID )						\
            __asm__ __volatile__ (									\
            "\n\t  movl $"#MARK_ID", %%ebx"	\
            "\n\t  .byte 0x64, 0x67, 0x90"	\
            : : : "memory" );

            #define IACA_UD_BYTES __asm__ __volatile__ ("\n\t .byte 0x0F, 0x0B");

        #else
            #define IACA_UD_BYTES { __asm _emit 0x0F \
                                    __asm _emit 0x0B}

            #define IACA_SSC_MARK(x) {__asm  mov ebx, x\
            __asm  _emit 0x64 \
            __asm  _emit 0x67 \
            __asm  _emit 0x90 }

            #define IACA_VC64_START __writegsbyte(111, 111);
            #define IACA_VC64_END   __writegsbyte(222, 222);

        #endif

        #define IACA_START {IACA_UD_BYTES \
                            IACA_SSC_MARK(111)}
        #define IACA_END {IACA_SSC_MARK(222) \
                          IACA_UD_BYTES}

    #endif

#endif