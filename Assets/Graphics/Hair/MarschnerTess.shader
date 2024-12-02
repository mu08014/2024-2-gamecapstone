// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Custom/MarschnerTess"
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
            #pragma target 3.0
            #pragma vertex vert
            #pragma geometry geom
            #pragma fragment frag
            #pragma hull hull
            #pragma domain domain
            #pragma require geometry
            #pragma require tessellation tessHW
            #pragma multi_compile_fwdbase
            #pragma multi_compile_shadowcaster
            #pragma addshadow
            #include "UnityCG.cginc"
            #include "AutoLight.cginc"

            #define MAX_HAIR_PARTICLE 15
            #define PARTICLE_COUNT 4
            #define BAZIER_COUNT 10
            #define _TessellationFactor 5
            #define _TessellationEdgeLength 1

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
                float2 uv : TEXCOORD1;

                float3 lightDir : TEXCOORD2;
                float3 cameraDir : TEXCOORD3;

                float3 middle : TEXCOORD5;
                float3 head : TEXCOORD6;
                float3 tail : TEXCOORD7;
                SHADOW_COORDS(8)
            };

            struct h2d_ConstantOutput {
                float TessEdge[3] : SV_TessFactor;
                float TessInside : SV_InsideTessFactor;
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                
                int hindex = (int)v.idx.x;
                int pindex = (int)v.idx.y;
                
                if (pindex == 0) {
                    o.head = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    o.tail = (_VertexPosition[hindex * PARTICLE_COUNT + pindex + 1] + o.head) / 2;
                    o.middle = (o.head + o.tail) / 2;
                }
                else if (pindex == PARTICLE_COUNT - 1) {
                    o.tail = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    o.head = (_VertexPosition[hindex * PARTICLE_COUNT + pindex - 1] + o.tail) / 2;
                    o.middle = (o.head + o.tail) / 2;
                } else {
                    o.middle = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    o.head = (_VertexPosition[hindex * PARTICLE_COUNT + pindex - 1] + o.middle) / 2;
                    o.tail = (_VertexPosition[hindex * PARTICLE_COUNT + pindex + 1] + o.middle) / 2;
                }
                
                o.uv = v.uv;
                o.pos = v.vertex;
                o.normal = v.normal;
                
                o.lightDir = -mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex);

                o.color = tex2Dlod(_MainTex, float4(v.uv, 0, 0));
                o.color = _Color;

                TRANSFER_SHADOW(o);
             
                return o;
            }

            // Hull Shader (Tessellation Control)
            [domain("tri")]
            [partitioning("integer")]
            [outputtopology("triangle_cw")]
            [outputcontrolpoints(3)]
            [patchconstantfunc("PatchConstantFunction")]
            //[maxtessfactor(64)]
            v2f hull(InputPatch<v2f, 3> input, 
                uint controlPointID : SV_OutputControlPointID)
            {
                return input[controlPointID];
            }

            h2d_ConstantOutput PatchConstantFunction(InputPatch<v2f, 3> patch) {
                h2d_ConstantOutput o = (h2d_ConstantOutput)0;
                float edge = _TessellationFactor * (1.0 / _TessellationEdgeLength);
                o.TessEdge[0] = o.TessEdge[1] = o.TessEdge[2] = edge;
                o.TessInside = edge;
                return o;
            }

            // Domain Shader (Tessellation Evaluation)
            [domain("tri")]
            v2f domain(
                h2d_ConstantOutput tessFactors, 
                OutputPatch<v2f, 3> patch,
                float3 bary : SV_DomainLocation) 
            {
                v2f output;
                
                output.head = bary.x * patch[0].head +
                            bary.y * patch[1].head +
                            bary.z * patch[2].head;
                             
                output.middle = bary.x * patch[0].middle +
                            bary.y * patch[1].middle +
                            bary.z * patch[2].middle;
                             
                output.tail = bary.x * patch[0].tail +
                            bary.y * patch[1].tail +
                            bary.z * patch[2].tail;

                output.normal = bary.x * patch[0].normal +
                            bary.y * patch[1].normal +
                            bary.z * patch[2].normal;

                output.lightDir = patch[0].lightDir;

                output.cameraDir = bary.x * patch[0].cameraDir +
                            bary.y * patch[1].cameraDir +
                            bary.z * patch[2].cameraDir;

                output.uv = bary.x * patch[0].uv +
                            bary.y * patch[1].uv +
                            bary.z * patch[2].uv;

                float3 pos = bary.x * patch[0].pos +
                             bary.y * patch[1].pos +
                             bary.z * patch[2].pos;

                output.pos = UnityObjectToClipPos(float4(pos, 1.0));
                return output;
            }
             
            // Geometry Shader: 메시 인덱스를 받아 삼각형으로 변환
            [maxvertexcount((BAZIER_COUNT + 1) * 2)]
            void geom(triangle v2f input[3], inout TriangleStream<v2f> triStream)
            {
                for (int i = 0; i < 1; i++) {
                    float3 pos0 = (input[0].head + input[1].head + input[2].head) / 3;
                    float3 pos1 = (input[0].middle + input[1].middle + input[2].middle) / 3;
                    float3 pos2 = (input[0].tail + input[1].tail + input[2].tail) / 3;

                    // set offset
                    float4 npos0 = UnityObjectToClipPos(float4(pos0, 1));
                    npos0 = npos0 / npos0.w;
                    float4 npos1 = UnityObjectToClipPos(float4(pos1, 1));
                    npos1 = npos1 / npos1.w;
                    float4 npos2 = UnityObjectToClipPos(float4(pos2, 1));
                    npos2 = npos2 / npos2.w;

                    float3 voffset0 = normalize(float3(npos0.xy - npos1.xy, 0));
                    float3 voffset1 = normalize(float3(npos1.xy - npos2.xy, 0));
                    float3 offset0 = cross(voffset0, float3(0, 0, 1));
                    float3 offset1 = cross(voffset1, float3(0, 0, 1));

                    // set tangent
                    float3 tan0 = normalize(pos0 - pos1);
                    float3 tan1 = normalize(pos1 - pos2);

                    for (int j = BAZIER_COUNT; j >= 0; j--) {
                        float3 n0 = (pos0 * j + pos1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        float3 n1 = (pos1 * j + pos2 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        float3 n = (n0 * j + n1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;

                        float4 nn = UnityObjectToClipPos(float4(n, 1));
                    
                        //float3 offset = float3(0.01, 0, 0);
                        float3 offset = (offset0 * j + offset1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        offset.y *= _AspectRatio;
                        offset *= 0.01;

                        v2f p0 = input[0];
                        p0.pos = nn;
                        p0.tangent = (tan0 * j + tan1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        p0.cameraDir = normalize(
                            (mul(unity_WorldToObject, _WorldSpaceCameraPos) - float4(n, 1)).xyz);

                        v2f v0 = p0;
                        v2f v1 = p0;
                        v0.pos.xy += offset.xy;
                        v1.pos.xy -= offset.xy;

                        TRANSFER_SHADOW(v0);
                        TRANSFER_SHADOW(v1);

                        triStream.Append(v0);
                        triStream.Append(v1);
                    }
                }
            }
            
            float4 frag (v2f i) : SV_Target
            {
                // diffuse
                float3 texColor = tex2D(_MainTex, i.uv);
                float DotTL = saturate(dot(i.tangent, i.lightDir));
                float3 diffuse = _DiffuseTerm * saturate(0.75 * sqrt(1 - DotTL * DotTL) + 0.25) * texColor;

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
                float3 embient = texColor * _AmbientTerm;
                
                float shadowAttenuation = SHADOW_ATTENUATION(i);
                float3 color = (diffuse + specular) * shadowAttenuation + embient;

                return float4(color, 1);
            }
            ENDCG
        }
        
        Pass
        {
            Cull Off
            Name "ShadowCaster"
            Tags { "LightMode" = "ShadowCaster" }

            CGPROGRAM
            #pragma target 3.0
            #pragma vertex vert
            #pragma geometry geom
            #pragma fragment frag
            #pragma hull hull
            #pragma domain domain
            #pragma require geometry
            #pragma require tessellation tessHW
            #pragma multi_compile_fwdbase
            #pragma multi_compile_shadowcaster
            #pragma addshadow
            #include "UnityCG.cginc"
            #include "AutoLight.cginc"

            #define MAX_HAIR_PARTICLE 15
            #define PARTICLE_COUNT 4
            #define BAZIER_COUNT 10
            #define _TessellationFactor 5
            #define _TessellationEdgeLength 1

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
                float2 uv : TEXCOORD1;

                float3 lightDir : TEXCOORD2;
                float3 cameraDir : TEXCOORD3;

                float3 middle : TEXCOORD5;
                float3 head : TEXCOORD6;
                float3 tail : TEXCOORD7;
                SHADOW_COORDS(8)
            };

            struct h2d_ConstantOutput {
                float TessEdge[3] : SV_TessFactor;
                float TessInside : SV_InsideTessFactor;
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                
                int hindex = (int)v.idx.x;
                int pindex = (int)v.idx.y;
                
                if (pindex == 0) {
                    o.head = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    o.tail = (_VertexPosition[hindex * PARTICLE_COUNT + pindex + 1] + o.head) / 2;
                    o.middle = (o.head + o.tail) / 2;
                }
                else if (pindex == PARTICLE_COUNT - 1) {
                    o.tail = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    o.head = (_VertexPosition[hindex * PARTICLE_COUNT + pindex - 1] + o.tail) / 2;
                    o.middle = (o.head + o.tail) / 2;
                } else {
                    o.middle = _VertexPosition[hindex * PARTICLE_COUNT + pindex];
                    o.head = (_VertexPosition[hindex * PARTICLE_COUNT + pindex - 1] + o.middle) / 2;
                    o.tail = (_VertexPosition[hindex * PARTICLE_COUNT + pindex + 1] + o.middle) / 2;
                }
                
                o.uv = v.uv;
                o.pos = v.vertex;
                o.normal = v.normal;
                
                o.lightDir = -mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex);

                o.color = tex2Dlod(_MainTex, float4(v.uv, 0, 0));
                o.color = _Color;

                TRANSFER_SHADOW(o);
             
                return o;
            }

            // Hull Shader (Tessellation Control)
            [domain("tri")]
            [partitioning("integer")]
            [outputtopology("triangle_cw")]
            [outputcontrolpoints(3)]
            [patchconstantfunc("PatchConstantFunction")]
            //[maxtessfactor(64)]
            v2f hull(InputPatch<v2f, 3> input, 
                uint controlPointID : SV_OutputControlPointID)
            {
                return input[controlPointID];
            }

            h2d_ConstantOutput PatchConstantFunction(InputPatch<v2f, 3> patch) {
                h2d_ConstantOutput o = (h2d_ConstantOutput)0;
                float edge = _TessellationFactor * (1.0 / _TessellationEdgeLength);
                o.TessEdge[0] = o.TessEdge[1] = o.TessEdge[2] = edge;
                o.TessInside = edge;
                return o;
            }

            // Domain Shader (Tessellation Evaluation)
            [domain("tri")]
            v2f domain(
                h2d_ConstantOutput tessFactors, 
                OutputPatch<v2f, 3> patch,
                float3 bary : SV_DomainLocation) 
            {
                v2f output;
                
                output.head = bary.x * patch[0].head +
                            bary.y * patch[1].head +
                            bary.z * patch[2].head;
                             
                output.middle = bary.x * patch[0].middle +
                            bary.y * patch[1].middle +
                            bary.z * patch[2].middle;
                             
                output.tail = bary.x * patch[0].tail +
                            bary.y * patch[1].tail +
                            bary.z * patch[2].tail;

                output.normal = bary.x * patch[0].normal +
                            bary.y * patch[1].normal +
                            bary.z * patch[2].normal;

                output.lightDir = patch[0].lightDir;

                output.cameraDir = bary.x * patch[0].cameraDir +
                            bary.y * patch[1].cameraDir +
                            bary.z * patch[2].cameraDir;

                output.uv = bary.x * patch[0].uv +
                            bary.y * patch[1].uv +
                            bary.z * patch[2].uv;

                float3 pos = bary.x * patch[0].pos +
                             bary.y * patch[1].pos +
                             bary.z * patch[2].pos;

                output.pos = UnityObjectToClipPos(float4(pos, 1.0));
                return output;
            }
             
            // Geometry Shader: 메시 인덱스를 받아 삼각형으로 변환
            [maxvertexcount((BAZIER_COUNT + 1) * 2)]
            void geom(triangle v2f input[3], inout TriangleStream<v2f> triStream)
            {
                for (int i = 0; i < 1; i++) {
                    float3 pos0 = (input[0].head + input[1].head + input[2].head) / 3;
                    float3 pos1 = (input[0].middle + input[1].middle + input[2].middle) / 3;
                    float3 pos2 = (input[0].tail + input[1].tail + input[2].tail) / 3;

                    // set offset
                    float4 npos0 = UnityObjectToClipPos(float4(pos0, 1));
                    npos0 = npos0 / npos0.w;
                    float4 npos1 = UnityObjectToClipPos(float4(pos1, 1));
                    npos1 = npos1 / npos1.w;
                    float4 npos2 = UnityObjectToClipPos(float4(pos2, 1));
                    npos2 = npos2 / npos2.w;

                    float3 voffset0 = normalize(float3(npos0.xy - npos1.xy, 0));
                    float3 voffset1 = normalize(float3(npos1.xy - npos2.xy, 0));
                    float3 offset0 = cross(voffset0, float3(0, 0, 1));
                    float3 offset1 = cross(voffset1, float3(0, 0, 1));

                    // set tangent
                    float3 tan0 = normalize(pos0 - pos1);
                    float3 tan1 = normalize(pos1 - pos2);

                    for (int j = BAZIER_COUNT; j >= 0; j--) {
                        float3 n0 = (pos0 * j + pos1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        float3 n1 = (pos1 * j + pos2 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        float3 n = (n0 * j + n1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;

                        float4 nn = UnityObjectToClipPos(float4(n, 1));
                    
                        //float3 offset = float3(0.01, 0, 0);
                        float3 offset = (offset0 * j + offset1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        offset.y *= _AspectRatio;
                        offset *= 0.01;

                        v2f p0 = input[0];
                        p0.pos = nn;
                        p0.tangent = (tan0 * j + tan1 * (BAZIER_COUNT - j)) / BAZIER_COUNT;
                        p0.cameraDir = normalize(
                            (mul(unity_WorldToObject, _WorldSpaceCameraPos) - float4(n, 1)).xyz);

                        v2f v0 = p0;
                        v2f v1 = p0;
                        v0.pos.xy += offset.xy;
                        v1.pos.xy -= offset.xy;

                        triStream.Append(v0);
                        triStream.Append(v1);
                    }
                }
            }
            
            float4 frag (v2f i) : SV_Target
            {
                return 0;
                /*
                // diffuse
                float3 texColor = tex2D(_MainTex, i.uv);
                float DotTL = saturate(dot(i.tangent, i.lightDir));
                float3 diffuse = _DiffuseTerm * saturate(0.75 * sqrt(1 - DotTL * DotTL) + 0.25) * texColor;

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
                float3 embient = i.color * _AmbientTerm * 100;
               
                TRANSFER_SHADOW(i);
                float shadowAttenuation = SHADOW_ATTENUATION(i);
                float3 color = (diffuse + specular) * shadowAttenuation  + embient;

                return float4(embient, 1);
                */
            }
            ENDCG
        }
        
        /*
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
        }*/
        
        //UsePass "VertexLit/SHADOWCASTER"
    }
    FallBack "Diffuse"
}
