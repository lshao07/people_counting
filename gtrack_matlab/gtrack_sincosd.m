function [pSinVal, pCosVal] = gtrack_sincosd(theta)
    pSinVal = sin(theta);
    pCosVal = cos(theta);
end

%{
void gtrack_sincosd(float theta, float * pSinVal, float * pCosVal)     
{     
    uint32_t i;                                /* Index for reading nearwst output values */     
    float x1 = -179.0f;                        /* Initial input value */     
    float ysin0, ysin1;                        /* nearest output values */     
    float ycos0, ycos1;                        /* nearest output values */     
    float fract;                               /* fractional part of input */     
    float ftheta = (float)((int32_t)theta);  


    /* index calculation for reading nearest output values */     
    i = (uint32_t) (theta - x1);     

    /* Calculation of fractional part */     
    if(theta > 0.0f)     
    {     
        fract = theta - ftheta;     
    }     
    else     
    {     
        fract = (theta - ftheta) + 1.0f;     
    }     

    /* reading nearest sine output values */     
    ysin0 = sinTable[i];     
    ysin1 = sinTable[i + 1u];     

    /* reading nearest cosine output values */     
    ycos0 = cosTable[i];     
    ycos1 = cosTable[i + 1u];     

    /* difference of nearest sine output value */  
    ysin1 = ysin1 - ysin0;  

    /* difference of nearest cosine output value */  
    ycos1 = ycos1 - ycos0;  

    /* Calculation of sine value */     
    *pSinVal = ysin0 + (fract * ysin1);

    /* Calculation of cosine value */     
    *pCosVal = ycos0 + (fract * ycos1);
}     
%}