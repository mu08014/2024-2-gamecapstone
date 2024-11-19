Shader "Custom/WaterRippleEffect"
{
    Properties
    {
        _MainTex ("Main Texture", 2D) = "white" {}
        _ColorTex ("Color Texture", 2D) = "white" {}
        _NoiseTex ("Noise Texture", 2D) = "white" {}
        _ReflectionTex ("Reflection Texture", 2D) = "white" {}
        _DistortionStrength ("Distortion Strength", Float) = 0.05
        _TimeScale ("Time Scale", Float) = 1.0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100
        
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            
            #include "UnityCG.cginc"
            
            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };
            
            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 pos : SV_POSITION;
            };
            
            sampler2D _MainTex;
            sampler2D _ColorTex;
            sampler2D _NoiseTex;
            sampler2D _ReflectionTex;
            float _DistortionStrength;
            float _TimeScale;
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                return o;
            }
            
            fixed4 frag (v2f i) : SV_Target
            {
                float2 noise = tex2D(_NoiseTex, i.uv).rg;
                if (noise.x * noise.x + noise.y + noise.y > 0.0001) {
                    noise *= 2.0f;
                    noise -= 1.0f;
                }
                float2 distortion = -noise * _DistortionStrength;
                float2 uv = i.uv + distortion;
                fixed4 target = tex2D(_MainTex, uv);

                fixed4 color = tex2D(_ColorTex, i.uv);

                // ssr
                float4 refl = tex2D(_ReflectionTex, i.uv);
                //float edgeFade = 1.0f - pow(length(refl.xy - 0.5f) * 2.0f, 1.0f);
                if (refl.w != 0) {
                    color.rgb = tex2D(_MainTex, refl.xy).rgb;
                    //color.a = 1;
                }
                return target * (1 - color.a) + color * color.a;
            }
            ENDCG
        }
    }
}
