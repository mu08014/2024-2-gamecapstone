// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Custom/KajiyaLine"
{
    Properties
    {
        _LineWidth("LineWidth", float) = 0.1
        _DiffuseTerm("DiffuseTerm", Range(0.0, 1.0)) = 1
        _SpecularTerm("SpecularTerm", Range(0.0, 1.0)) = 1
        _SpecularPower("SpecularPower", Range(0.0, 1000.0)) = 10
        _Color ("Color", Color) = (1,1,1,1)
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
            #pragma fragment frag
            #pragma multi_compile_fwdbase
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
                float3 color : COLOR;
                float4 pos : SV_POSITION;
                float3 tangent : TANGENT;
                float3 normal : NORMAL;

                float4 lightDir : TEXCOORD3;
                float4 cameraDir : TEXCOORD2;

                SHADOW_COORDS(1)
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.normal = v.normal;
                o.tangent = normalize(v.tangent);
                
                o.lightDir = mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex);

                o.color = tex2Dlod(_MainTex, float4(v.uv, 0, 0));

                TRANSFER_SHADOW(o);
                
                return o;
            }

            float4 frag (v2f i) : SV_Target
            {
                float diffuse_val = saturate(dot(i.lightDir, i.normal));
                float3 diffuse = i.color;

                float3 halfVec = normalize(i.lightDir + i.cameraDir);
                float3 specular = pow( 
                    sqrt(1 - dot(i.tangent, halfVec) * dot(i.tangent, halfVec)),
                    _SpecularPower) * _SpecularTerm;
                  
                float3 embient = i.color * 0.2f;

                float3 color = (diffuse + specular) * diffuse_val  + embient;

                float shadowAttenuation = SHADOW_ATTENUATION(i);
                color = color * shadowAttenuation;

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

                // ShadowCaster Pass 추가
        Pass
        {
            Name "ShadowCaster"
            Tags { "LightMode" = "ShadowCaster" }

            CGPROGRAM
            #pragma vertex vert
            #pragma geometry geom
            #pragma fragment frag

            // 그림자 렌더링용 쉐이더 코드를 포함
            #pragma multi_compile_shadowcaster

            struct v2f
            {
                float4 pos : SV_POSITION;
            };

            // 기본 Vertex Shader (ShadowCaster 패스)
            v2f vert(appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                return o;
            }

            // Geometry Shader (ShadowCaster 패스에서도 Geometry Shader 사용)
            [maxvertexcount(6)]
            void geom(line v2f input[2], inout TriangleStream<v2f> triStream)
            {
                float3 cameraRight = normalize(cross(_WorldSpaceCameraPos.xyz - input[0].pos.xyz, float3(0, 1, 0)));
                float3 offset = cameraRight * _LineWidth;

                // 첫 번째 삼각형
                v2f v0 = input[0];
                v2f v1 = input[1];
                v2f v2 = input[0];
                v0.pos.xy += offset.xy;
                v1.pos.xy += offset.xy;
                v2.pos.xy -= offset.xy;

                triStream.Append(v0);
                triStream.Append(v1);
                triStream.Append(v2);

                // 두 번째 삼각형
                v2f v3 = input[1];
                v2f v4 = input[1];
                v2f v5 = input[0];
                v3.pos.xy -= offset.xy;
                v4.pos.xy -= offset.xy;
                v5.pos.xy += offset.xy;

                triStream.Append(v3);
                triStream.Append(v4);
                triStream.Append(v5);
            }

            half4 frag(v2f i) : SV_Target
            {
                return float4(0,0,0,1); // 그림자용 검정색 출력
            }
            ENDCG
        }

        UsePass "VertexLit/SHADOWCASTER"
    }
    FallBack "Diffuse"
}
