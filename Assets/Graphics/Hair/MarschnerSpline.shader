// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Custom/MarschnerSpline"
{
    Properties
    {
        _AspectRatio("AspectRatio", float) = 1.0
        _LineWidth("LineWidth", float) = 0.1
        _DiffuseTerm("DiffuseTerm", Range(0.0, 1.0)) = 1
        _SpecularRTerm("SpecularRTerm", Range(0.0, 1.0)) = 1
        _SpecularTTTerm("SpecularTTTerm", Range(0.0, 1.0)) = 1
        _SpecularTRTTerm("SpecularTRTTerm", Range(0.0, 1.0)) = 1
        _AmbientTerm("AmbientTerm", Range(0.0, 1.0)) = 1
        _SpecularPower("SpecularPower", Range(0.0, 1000.0)) = 10
        _Color ("Color", Color) = (1,1,1,1)
        _SecondSpecularColor("SecondColor", Color) = (1,1,1,1)
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" "Queue"="Geometry" }
        
        Pass
        {
            Cull Off

            CGPROGRAM
            //#pragma target 3.0
            #pragma vertex vert
            #pragma geometry geom
            #pragma fragment frag
            #pragma multi_compile_fwdbase
            #pragma multi_compile_shadowcaster
            #include "UnityCG.cginc"
            #include "AutoLight.cginc"

            #define MAX_HAIR_PARTICLE 15
            #define PARTICLE_COUNT 4
            #define BAZIER_COUNT 19

            StructuredBuffer<float3> _VertexPosition; // 모든 position을 한 번에 전달
            
            float _AspectRatio;
            float _LineWidth;
            float _DiffuseTerm;
            float _SpecularRTerm; 
            float _SpecularTTTerm;
            float _SpecularTRTTerm;
            float _AmbientTerm;
            float _SpecularPower;
            fixed4 _Color;
            fixed4 _SecondSpecularColor;
            sampler2D _MainTex;

            struct appdata
            {
                float4 vertex : POSITION;
                float3 tangent : TEXCOORD1;
                float3 normal : TEXCOORD2;
                float2 uv : TEXCOORD3;
                float2 idx: TEXCOORD4;
            };
            
            struct v2f
            {
                float3 color : COLOR;
                float4 pos : SV_POSITION;
                float3 tangent : TANGENT;
                float3 normal : NORMAL;

                float3 lightDir : TEXCOORD2;
                float3 cameraDir : TEXCOORD3;
                float2 idx : TEXCOORD4;

                SHADOW_COORDS(1)
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = float4(0, 0, 0, 0);
                o.normal = v.normal;
                o.tangent = normalize(v.tangent);
                
                o.lightDir = -mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex);

                o.color = tex2Dlod(_MainTex, float4(v.uv, 0, 0));
                o.color = _Color;

                o.idx = v.idx;

                TRANSFER_SHADOW(o);
             
                return o;
            }
             
             
            // Geometry Shader: 메시 인덱스를 받아 삼각형으로 변환
            [maxvertexcount((BAZIER_COUNT + 1) * 2)]
            void geom(point v2f input[1], inout TriangleStream<v2f> triStream)
            { 
                int hindex = (int)input[0].idx.x;
                int pindex = (int)input[0].idx.y;

                float3 pos0;
                float3 pos1;
                float3 pos2;

                if (pindex == 0) {
                    pos0 = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    pos2 = (_VertexPosition[hindex * PARTICLE_COUNT + pindex + 1] + pos0) / 2.0;
                    pos1 = (pos0 + pos2) / 2.0;
                } else if (pindex == PARTICLE_COUNT - 1) {
                    pos0 = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    pos2 = (_VertexPosition[hindex * PARTICLE_COUNT + pindex - 1] + pos0) / 2.0;
                    pos1 = (pos0 + pos2) / 2.0;
                } else {
                    pos1 = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    pos0 = (_VertexPosition[hindex * PARTICLE_COUNT + pindex - 1] + pos1) / 2.0;
                    pos2 = (_VertexPosition[hindex * PARTICLE_COUNT + pindex + 1] + pos1) / 2.0;
                }
                    
                for (int j = BAZIER_COUNT; j >= 0; j--) {
                    float3 n0 = (pos0 * j + pos1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                    float3 n1 = (pos1 * j + pos2 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                    float3 n = (n0 * j + n1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;

                    float4 nn = UnityObjectToClipPos(float4(n, 1));
                    
                    float3 offset = float3(0.01, 0, 0);
                    v2f p0 = input[0]; p0.pos = nn;
                    v2f v0 = p0, v1 = p0;
                    v0.pos.xy += offset.xy;
                    v1.pos.xy -= offset.xy;

                    triStream.Append(v0);
                    triStream.Append(v1);
                }
            }

            float4 frag (v2f i) : SV_Target
            {
                return float4(1, 0, 0, 1);

                // diffuse
                float DotTL = saturate(dot(i.tangent, i.lightDir));
                float3 diffuse = _DiffuseTerm * saturate(0.75 * sqrt(1 - DotTL * DotTL) + 0.25) * i.color;

                float DotTC = max(dot(i.tangent, i.cameraDir), 0.0);
                float SinTL = sqrt(1 - DotTL * DotTL);
                float SinTC = sqrt(1 - DotTC * DotTC);

                float a = 3.0 * 3.14159265359 / 180.0; 

                // specular
                float3 specularR = _SpecularRTerm * pow(
                //    DotTL * DotTC + SinTL * SinTC,
                    (DotTL*cos(2*a) - SinTL*sin(2*a)) * DotTC + (SinTL*cos(2*a) + DotTL*sin(2*a)) * SinTC,
                    _SpecularPower / 2);

                    
                float3 specularTRT = _SpecularTRTTerm * pow(
                    (DotTL*cos(3*a) + SinTL*sin(3*a)) * DotTC + (SinTL*cos(3*a) - DotTL*sin(3*a)) * SinTC,
                    _SpecularPower) * _SecondSpecularColor;

                float3 specular = specularR + specularTRT;
                
                // embient
                float3 embient = i.color * _AmbientTerm;
                
                float shadowAttenuation = SHADOW_ATTENUATION(i);
                float3 color = (diffuse + specular) * shadowAttenuation  + embient;

                //return float4(h.xyz, 1);
                return float4(color, 1);
            }
            ENDCG
        }

        Pass
        {
            Tags { "LightMode" = "ShadowCaster" }

            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_shadowcaster
            #include "UnityCG.cginc"
            #include "AutoLight.cginc"

            float _LineWidth;
            float _DiffuseTerm;
            float _SpecularTerm;
            float _SpecularPower;
            fixed4 _Color;
            sampler2D _MainTex;
            
            struct appdata
            {
                float4 vertex : POSITION;
                float3 tangent : TEXCOORD1;
                float3 normal : TEXCOORD2;
                float2 uv : TEXCOORD3;
            };

            struct v2f
            {
                V2F_SHADOW_CASTER;
            };

            v2f vert(appdata_base v)
            {
                v2f o;
                TRANSFER_SHADOW_CASTER(o)
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                return 0;
            }
            ENDCG
        }
        //UsePass "VertexLit/SHADOWCASTER"
    }
    FallBack "Diffuse"
}
